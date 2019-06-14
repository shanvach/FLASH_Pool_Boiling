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

  use Grid_interface, only  : Grid_getBlkPtr, Grid_releaseBlkPtr

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  integer, intent(in) :: blockID
  !!$ ---------------------------------
 
  real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData,facezData

  ! Point to Blocks centered variables:
  call Grid_getBlkPtr(blockID,solnData,CENTER)

  ! Point to Blocks face variables: 
  call Grid_getBlkPtr(blockID,facexData,FACEX)
  call Grid_getBlkPtr(blockID,faceyData,FACEY)

#if NDIM == 3
  call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

  solnData(PRES_VAR,:,:,:) = 0.0
  solnData(DELP_VAR,:,:,:) = 0.0
  solnData(DUST_VAR,:,:,:) = 0.0
  solnData(TVIS_VAR,:,:,:) = 0.0
  solnData(TEMP_VAR,:,:,:) = 0.5

  facexData(VELC_FACE_VAR,:,:,:) = 1.0
  faceyData(VELC_FACE_VAR,:,:,:) = 0.0
  facexData(RHDS_FACE_VAR,:,:,:) = 0.0
  faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
  facezData(VELC_FACE_VAR,:,:,:) = 0.0
  facezData(RHDS_FACE_VAR,:,:,:) = 0.0
#endif

  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

#if NDIM == 3
  call Grid_releaseBlkPtr(blockID,solnData,FACEZ)
#endif

  return

111    format (i4,3x,i4)
112    format (3(3x,e12.4))

end subroutine Simulation_initBlock
