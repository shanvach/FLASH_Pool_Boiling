subroutine sim_trackBubble(blockCount, blockList,vol,x,y,vx,vy,xh,yh,xl,yl)

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

  use Driver_data, ONLY : dr_simTime, dr_meshMe

  implicit none

#include "constants.h"
#include "Flash.h"
include "Flash_mpi.h"

  !!$ Arguments -----------------------
  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList !blockCount
  real :: vol, x, y, vx, vy, xh, yh, xl, yl 
  !!$ ---------------------------------
 
  integer :: i, j, k
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) ::  blIndSize,blIndSizeGC

  real, dimension(MDIM)  :: coord,bsize,del
  real ::  boundBox(2,MDIM)
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData

  real :: xcell,xedge,ycell,yedge
  real :: volx, voly, volvx, volvy, u, v
  real :: volAll, volxAll, volyAll, volvxAll, volvyAll
  real :: xhl, yhl, xll, yll, yf

  integer :: lb, blockID, ierr

  
  vol = 0.0d0
  volx = 0.0d0
  voly = 0.0d0
  volvx = 0.0d0
  volvy = 0.0d0
  
  yhl = -1.0d30
  yll = 1.0d30
  xhl = -1.0d30
  xll = -1.0d30

  do lb = 1,blockCount
    blockID = blockList(lb)

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

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)

!  write(1000000+1000*dr_meshMe+lb,*) 'coord', coord
!  write(1000000+1000*dr_meshMe+lb,*) 'bsize', bsize
!  write(1000000+1000*dr_meshMe+lb,*) 'del', del
!  write(1000000+1000*dr_meshMe+lb,*) 'blkLimits', blkLimits, blkLimitsGC
 
  !- kpd - Initialize the distance function in the 1st quadrant 
  !do k=1,blkLimitsGC(HIGH,KAXIS)
     do j=blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
           if(solnData(DFUN_VAR,i,j,1) > 0.0d0) then
           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS)

           ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS)  +  &
                   0.5*del(JAXIS)

           u = 0.5d0*(facexData(VELC_FACE_VAR,i,j,1) + facexData(VELC_FACE_VAR,i+1,j,1))
           v = 0.5d0*(faceyData(VELC_FACE_VAR,i,j,1) + faceyData(VELC_FACE_VAR,i,j+1,1))
           vol = vol + del(IAXIS)*del(JAXIS)
           volx = volx + xcell*del(IAXIS)*del(JAXIS)
           voly = voly + ycell*del(IAXIS)*del(JAXIS)
           volvx = volvx + u*del(IAXIS)*del(JAXIS)
           volvy = volvy + v*del(IAXIS)*del(JAXIS)          
!           write(1000000+1000*dr_meshMe+lb,*) 'DFUN', solnData(DFUN_VAR,i,j,1), i, j, vol
           endif
       
           ! Only for free surface flow, in which the gas phase is above the liquid phase
           if(solnData(DFUN_VAR,i,j,1) > 0.0d0 .and. solnData(DFUN_VAR,i,j-1,1) < 0.0d0) then
             xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS)

             ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS)  +  &
                   0.5*del(JAXIS)
             yf = ycell - solnData(DFUN_VAR,i,j,1)

             if(yf < yll) yll = yf
             if(yf > yhl) yhl = yf
           endif
           !if(abs(solnData(DFUN_VAR,i,j,1)) <= sqrt(del(IAXIS)**2+del(JAXIS)**2)/2.0d0) then
           !xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
           !        real(i - NGUARD - 1)*del(IAXIS) +   &
           !        0.5*del(IAXIS)

           !ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
           !        real(j - NGUARD - 1)*del(JAXIS)  +  &
           !        0.5*del(JAXIS)

           !  if(ycell < yll) then
           !    yll = ycell
           !    write(1000000+1000*dr_meshMe+lb,*) 'yll', yll, xcell, solnData(DFUN_VAR,i,j,1)
           !  endif
           !  if(ycell > yhl) then
           !    yhl = ycell
           !    write(2000000+1000*dr_meshMe+lb,*) 'yhl', yhl, xcell, solnData(DFUN_VAR,i,j,1)
           !  endif
  !           if(ycell > yhl) yhl = ycell
           !endif
        enddo
     enddo
!  enddo

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

  enddo


  call MPI_Allreduce(yhl, yh, 1, FLASH_REAL,&
                     MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(yll, yl, 1, FLASH_REAL,&
                     MPI_MIN, MPI_COMM_WORLD, ierr)

! =====

  do lb = 1,blockCount
    blockID = blockList(lb)

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

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)

!  write(1000000+1000*dr_meshMe+lb,*) 'coord', coord
!  write(1000000+1000*dr_meshMe+lb,*) 'bsize', bsize
!  write(1000000+1000*dr_meshMe+lb,*) 'del', del
!  write(1000000+1000*dr_meshMe+lb,*) 'blkLimits', blkLimits, blkLimitsGC
 
  !- kpd - Initialize the distance function in the 1st quadrant 
  !do k=1,blkLimitsGC(HIGH,KAXIS)
     do j=blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
           if(solnData(DFUN_VAR,i,j,1) > 0.0d0 .and. solnData(DFUN_VAR,i,j-1,1) < 0.0d0) then
             xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS)

             ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS)  +  &
                   0.5*del(JAXIS)
             yf = ycell - solnData(DFUN_VAR,i,j,1)

             if(abs(yf - yll)<1.0d-4) xll = xcell
             if(abs(yf - yhl)<1.0d-4) xhl = xcell
           endif
           !if(abs(solnData(DFUN_VAR,i,j,1)) <= sqrt(del(IAXIS)**2+del(JAXIS)**2)/2.0d0) then
           !xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
           !        real(i - NGUARD - 1)*del(IAXIS) +   &
           !        0.5*del(IAXIS)

           !ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
            !       real(j - NGUARD - 1)*del(JAXIS)  +  &
            !       0.5*del(JAXIS)

            ! if(abs(ycell - yl)<del(JAXIS)) then
            !     xll = xcell
            !     write(3000000+1000*dr_meshMe+lb,*) 'xll', xll, yl, solnData(DFUN_VAR,i,j,1)
            ! endif
            ! if(abs(ycell - yh)<del(JAXIS)) then
            !     xhl = xcell
            !     write(4000000+1000*dr_meshMe+lb,*) 'xhl', xhl, yh, solnData(DFUN_VAR,i,j,1)
            ! endif
             !if(abs(ycell - yh)<del(JAXIS)) xhl = xcell
           !endif
        enddo
     enddo
!  enddo

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

  enddo

  call MPI_Allreduce(vol, volAll, 1, FLASH_REAL,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(volx, volxAll, 1, FLASH_REAL,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(voly, volyAll, 1, FLASH_REAL,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(volvx, volvxAll, 1, FLASH_REAL,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(volvy, volvyAll, 1, FLASH_REAL,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)

  call MPI_Allreduce(xhl, xh, 1, FLASH_REAL,&
                     MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(xll, xl, 1, FLASH_REAL,&
                     MPI_MAX, MPI_COMM_WORLD, ierr)

!  write(1000000+1000*dr_meshMe+lb,*) 'volAll', volAll

  vol = volAll
  x = volxAll/volAll
  y = volyAll/volAll
  vx = volvxAll/volAll
  vy = volvyAll/volAll

  return

end subroutine sim_trackBubble
