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

  real :: xcell,xedge,ycell,yedge
  real :: u, v, p

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


  facexData(VELC_FACE_VAR,:,:,:) = 1.0
  faceyData(VELC_FACE_VAR,:,:,:) = 0.0
  facexData(RHDS_FACE_VAR,:,:,:) = 0.0
  faceyData(RHDS_FACE_VAR,:,:,:) = 0.0


  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)

  ! Initial solution for U velocities:
   do j=1,blkLimitsGC(HIGH,JAXIS)
     do i=1,blkLimitsGC(HIGH,IAXIS)+1
       
       xedge  = coord(IAXIS) - bsize(IAXIS)/2.0 +  &
                 real(i - NGUARD - 1)*del(IAXIS)  
       ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                 real(j - NGUARD - 1)*del(JAXIS)  +  &
                 0.5*del(JAXIS)
       call sm_setPotFlowEllp(xedge, ycell,u,v,p) 
       facexData(VELC_FACE_VAR,i,j,1) = u
     enddo
   enddo
  ! Initial solution for V velocities:
   do j=1,blkLimitsGC(HIGH,JAXIS) +1
     do i=1,blkLimitsGC(HIGH,IAXIS)
       
       xcell  = coord(IAXIS) - bsize(IAXIS)/2.0 +  &
                 real(i - NGUARD - 1)*del(IAXIS) + &
                0.5*del(IAXIS)  
       yedge  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                 real(j - NGUARD - 1)*del(JAXIS)  
       call sm_setPotFlowEllp(xcell, yedge,u,v,p) 
       faceyData(VELC_FACE_VAR,i,j,1) = v
     enddo
   enddo
  ! Initial solution for pressure:
   do j=1,blkLimitsGC(HIGH,JAXIS)
     do i=1,blkLimitsGC(HIGH,IAXIS)
       
       xcell  = coord(IAXIS) - bsize(IAXIS)/2.0 +  &
                 real(i - NGUARD - 1)*del(IAXIS) + &
                0.5*del(IAXIS)  
       ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                 real(j - NGUARD - 1)*del(JAXIS) + &
                0.5*del(JAXIS) 
       call sm_setPotFlowEllp(xcell, ycell,u,v,p) 
       solnData(PRES_VAR,i,j,1) = p
     enddo
   enddo

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)



  return

111    format (i4,3x,i4)
112    format (3(3x,e12.4))

end subroutine Simulation_initBlock
