!!****if* source/Simulation/SimulationMain/SolidMechanics/First_test_Rigid/Simulation_initBlock
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


  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkBoundBox,    &
                             Grid_getBlkCenterCoords
  implicit none
#include "Flash.h"
#include "constants.h"

  !!$ Arguments -----------------------
  integer, intent(in) :: blockID
  !!$ ---------------------------------

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real ::  boundBox(2,MDIM)

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  integer :: i,j,k
  real :: coord(MDIM),bsize(MDIM),del(MDIM),xcell,ycell,zcell,xedge,yedge,zedge
  real :: xp_u(MDIM),xp_v(MDIM),xp_w(MDIM)

  ! Get Coord and Bsize for the block:
  ! Bounding box:
  call Grid_getBlkBoundBox(blockId,boundBox)
  bsize(:) = boundBox(2,:) - boundBox(1,:)

  call Grid_getBlkCenterCoords(blockId,coord)

  ! Get blocks dx, dy ,dz:
  call Grid_getDeltas(blockID,del)

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)

  ! Point to Blocks centered variables:
  call Grid_getBlkPtr(blockID,solnData,CENTER)

  ! Point to Blocks face variables: 
  call Grid_getBlkPtr(blockID,facexData,FACEX)
  call Grid_getBlkPtr(blockID,faceyData,FACEY)
  call Grid_getBlkPtr(blockID,facezData,FACEZ)

  solnData(PRES_VAR,:,:,:) = 1.

  ! U velocity:
  do k=1,blkLimitsGC(HIGH,KAXIS)

     zedge  = coord(KAXIS) - bsize(KAXIS)/2.0  +  &
              real(k - NGUARD - 1)*del(KAXIS)  

     zcell  = zedge + 0.5*del(KAXIS)

     do j=1,blkLimitsGC(HIGH,JAXIS)

        yedge  = coord(JAXIS) - bsize(JAXIS)/2.0  +  &
                 real(j - NGUARD - 1)*del(JAXIS) 

        ycell  = yedge + 0.5*del(JAXIS)

        do i=1,blkLimitsGC(HIGH,IAXIS)

           xedge  = coord(IAXIS) - bsize(IAXIS)/2.0  +  &
                    real(i - NGUARD - 1)*del(IAXIS)

           xcell  = xedge + 0.5*del(IAXIS)

           xp_u(1:MDIM) = (/ xedge, ycell, zcell /)
           xp_v(1:MDIM) = (/ xcell, yedge, zcell /)
           xp_w(1:MDIM) = (/ xcell, ycell, zedge /)

           facexData(VELC_FACE_VAR,i,j,k) = -1./sqrt(3.)*sum(xp_u(1:MDIM))* &
                                                         0.5/sqrt(1.5)
           faceyData(VELC_FACE_VAR,i,j,k) = -1./sqrt(3.)*sum(xp_v(1:MDIM))* & 
                                                         0.5/sqrt(1.5)
           facezData(VELC_FACE_VAR,i,j,k) =  1./sqrt(3.)*sum(xp_w(1:MDIM))* &
                                                         1.0/sqrt(1.5)

        enddo
     enddo
  enddo

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
  call Grid_releaseBlkPtr(blockID,facezData,FACEZ)


 
  return

end subroutine Simulation_initBlock
