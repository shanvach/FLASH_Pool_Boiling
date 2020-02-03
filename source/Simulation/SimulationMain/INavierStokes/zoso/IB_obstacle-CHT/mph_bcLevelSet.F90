! Set INFLOW and OUTFLOW boundary conditions for level set function
! Shizhao Wang
! Jun 12, 2015

!subroutine mph_bcLevelSet(blockID,u,t)
subroutine mph_bcLevelSet(u,t)

  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkBoundBox,    &
                             Grid_getBlkCenterCoords, &
                             Grid_getBlkBC

  use Driver_data, ONLY : dr_simTime, dr_meshMe

implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  real, intent(in) :: u, t
  !!$ ---------------------------------
 
  integer :: i, j, k, lb, count
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC, blkBC

  integer ,dimension(MAXBLOCKS) :: blkList
  integer :: blockID

  real, dimension(MDIM)  :: coord,bsize,del
  real ::  boundBox(2,MDIM)
  real, pointer, dimension(:,:,:,:) :: solnData

  real :: xcell, xedge, ycell, yedge

  call Grid_getListOfBlocks(ALL_BLKS, blkList,count)

!  write(200000+dr_meshMe,*) 'No. blks on', dr_meshMe, ':', count
!  write(200000+dr_meshMe,*) 'list', blkList(1:count)
!  write(200000+dr_meshMe,*) '===='

do lb = 1, count
  blockID = blkList(lb)

  ! Get Coord and Bsize for the block:
  ! Bounding box:
  call Grid_getBlkBoundBox(blockId,boundBox)
  bsize(:) = boundBox(2,:) - boundBox(1,:)

  call Grid_getBlkCenterCoords(blockId,coord)

  ! Get blocks dx, dy ,dz:
  call Grid_getDeltas(blockID,del)

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)

  call Grid_getBlkBC(blockID, blkBC)
!  write(300000+dr_meshMe*1000+lb,*) 'blk BC', blkBC(:,IAXIS), 'pos', coord

  ! Point to Blocks centered variables:
  call Grid_getBlkPtr(blockID,solnData,CENTER)

  ! Set the level set founction on the boundary 
   if(blkBC(LOW,IAXIS) == INFLOW_INS) then
!     write(400000+dr_meshMe*1000+lb,*) "VARIABLES = 'x', 'y', 'phi'"
!     write(400000+dr_meshMe*1000+lb,*) "ZONE I =", blkLimits(LOW,IAXIS)-1, ',J =', blkLimitsGC(HIGH,JAXIS) 
     k = 1
     do j=1,blkLimitsGC(HIGH,JAXIS)
        do i=1,blkLimits(LOW,IAXIS)-1

           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS) - u*t

           ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS)  +  &
                   0.5*del(JAXIS)

          solnData(DFUN_VAR,i,j,k) = ycell
 !         solnData(DFUN_VAR,i,j,k) = ycell-SIN(xcell)
           
 !         write(400000+dr_meshMe*1000+lb,*) xcell, ycell, solnData(DFUN_VAR,i,j,k)
        enddo
     enddo
   endif

   if(blkBC(HIGH,IAXIS) == OUTFLOW_INS) then
!     write(410000+dr_meshMe*1000+lb,*) "VARIABLES = 'x', 'y', 'phi'"
!     write(410000+dr_meshMe*1000+lb,*) "ZONE I =", blkLimits(LOW,IAXIS)-1, ',J =', blkLimitsGC(HIGH,JAXIS) 
     k = 1
     do j=1,blkLimitsGC(HIGH,JAXIS)
        do i=blkLimits(HIGH,IAXIS)+1, blkLimitsGC(HIGH,IAXIS)

           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS) - u*t

           ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS)  +  &
                   0.5*del(JAXIS)

          solnData(DFUN_VAR,i,j,k) = ycell
          !solnData(DFUN_VAR,i,j,k) = ycell-SIN(xcell)
          
!          write(410000+dr_meshMe*1000+lb,*) xcell, ycell, solnData(DFUN_VAR,i,j,k)

        enddo
     enddo
   endif

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

enddo

  return
endSubroutine mph_bcLevelSet

