!!****if* source/Simulation/SimulationMain/INavierStokes/3D/ChannelLam/Simulation_initBlock
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

  use Simulation_data, ONLY : sim_xMin, sim_xMax, &
                              sim_yMin, sim_yMax, &
                              sim_zMin, sim_zMax, &
                              sim_gCell

  use Grid_interface, ONLY : Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  integer, intent(in) :: blockID
  !!$ ---------------------------------
 
  integer :: i, j, k
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) ::  blIndSize,blIndSizeGC

  real, dimension(MDIM)  :: coord,bsize
  real ::  boundBox(2,MDIM)
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  !----------------------------------------------------------------------


  ! Point to Blocks centered variables:
  call Grid_getBlkPtr(blockID,solnData,CENTER)

  ! Point to Blocks face variables: 
  call Grid_getBlkPtr(blockID,facexData,FACEX)
  call Grid_getBlkPtr(blockID,faceyData,FACEY)
  call Grid_getBlkPtr(blockID,facezData,FACEZ)

!  CALL RANDOM_SEED

  ! set values for u,v velocities and pressure
  ! In this problem we start with 0 velocities and leave a the unit body force in
  ! the x direction to drive the flow.

!  call random_number(facexData(VELC_FACE_VAR,:,:,:))
!  call random_number(faceyData(VELC_FACE_VAR,:,:,:))
!  call random_number(facezData(VELC_FACE_VAR,:,:,:))


  facexData(VELC_FACE_VAR,:,:,:) = 1.0
  faceyData(VELC_FACE_VAR,:,:,:) = 0.01
  facezData(VELC_FACE_VAR,:,:,:) = 0.01
  solnData(PRES_VAR,:,:,:) = 0.0
  solnData(DELP_VAR,:,:,:) = 0.0
  solnData(DUST_VAR,:,:,:) = 0.0
  solnData(TVIS_VAR,:,:,:) = 0.0

  facexData(RHDS_FACE_VAR,:,:,:) = 0.0
  faceyData(RHDS_FACE_VAR,:,:,:) = 0.0
  facezData(RHDS_FACE_VAR,:,:,:) = 0.0


  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
  call Grid_releaseBlkPtr(blockID,facezData,FACEZ)


  return

end subroutine Simulation_initBlock
