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

  use Simulation_data, ONLY : sim_xMin, sim_xMax, &
                              sim_yMin, sim_yMax, &
                              sim_gCell

  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkBoundBox,    &
                             Grid_getBlkCenterCoords

  use Driver_data, ONLY : dr_simTime

  ! for all PERIODIC boundary conditions, Shizhao Wang, Nov 2014
  use SolidMechanics_data, ONLY : sm_gravX,sm_gravY,sm_gravZ

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  integer, intent(in) :: blockID
  !!$ ---------------------------------
 
  integer :: i, j, k
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) ::  blIndSize,blIndSizeGC

  real, dimension(MDIM)  :: coord,bsize,del
  real ::  boundBox(2,MDIM)
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData

  real :: xcell,xedge,ycell,yedge,zcell, zedge


  !----------------------------------------------------------------------
  
  ! Get Coord and Bsize for the block:
  ! Bounding box:
  call Grid_getBlkBoundBox(blockId,boundBox)
  bsize(:) = boundBox(2,:) - boundBox(1,:)

  call Grid_getBlkCenterCoords(blockId,coord)

  ! Get blocks dx, dy ,dz:
  call Grid_getDeltas(blockID,del)

  ! Point to Blocks centered variables:
  call Grid_getBlkPtr(blockID,solnData,CENTER)

  ! Point to Blocks face variables: 
  call Grid_getBlkPtr(blockID,facexData,FACEX)
  call Grid_getBlkPtr(blockID,faceyData,FACEY)


  ! set values for u,v velocities and pressure
  solnData(PRES_VAR,:,:,:) = 0.0
  solnData(DELP_VAR,:,:,:) = 0.0
  solnData(DUST_VAR,:,:,:) = 0.0
  solnData(TVIS_VAR,:,:,:) = 0.0


  facexData(VELC_FACE_VAR,:,:,:) = 0.0
  faceyData(VELC_FACE_VAR,:,:,:) = 0.0
  facexData(RHDS_FACE_VAR,:,:,:) = 0.0
  faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

  ! For the flows with all PERIODIC boundary conditions
  !  Shizhao Wang, Nov 2014
  !call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  !write(*,*) 'gravity:', sm_gravX, sm_gravY, sm_gravZ, 'in initializing block', blockID
  !write(*,*) 'blockIndLimit:', blkLimits, 'in initializing block', blockID

!#if NDIM == MDIM
!  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
!    zcell = coord(KAXIS) - bsize(KAXIS)/2 + del(KAXIS)/2 + (k-NGUARD-1)*del(KAXIS)
!    do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
!      ycell = coord(JAXIS) - bsize(JAXIS)/2 + del(JAXIS)/2 + (j-NGUARD-1)*del(JAXIS)
!      do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
!        xcell = coord(IAXIS) - bsize(IAXIS)/2 + del(IAXIS)/2 + (i-NGUARD-1)*del(IAXIS)
!        solnData(PRES_VAR,i,j,k) = sm_gravX*xcell + sm_gravY*ycell + sm_gravZ*zcell
!      enddo
!    enddo
!  enddo
!#else
!  do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS) 
!    ycell = coord(JAXIS) - bsize(JAXIS)/2 + del(JAXIS)/2 + (j-NGUARD-1)*del(JAXIS)
!    do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
!      xcell = coord(IAXIS) - bsize(IAXIS)/2 + del(IAXIS)/2 + (i-NGUARD-1)*del(IAXIS)
!      solnData(PRES_VAR,i,j,1) = sm_gravX*xcell + sm_gravY*ycell
!    enddo
!  enddo
!#endif
 
  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)



  return

111    format (i4,3x,i4)
112    format (3(3x,e12.4))

end subroutine Simulation_initBlock
