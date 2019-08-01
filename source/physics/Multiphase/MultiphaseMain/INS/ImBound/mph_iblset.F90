!! source/physics/Multiphase/MultiphaseMain/INS/ImBound
!!
!! NAME
!!
!! mph_iblset(blockCount,blockList)
!!
!! SYNOPSIS
!! 
!!  
!! VARIABLES
!!
!!
!! DESCRIPTION
!! 
!! Subroutine to find the distance function lambda for 
!! the immersed boundary (IB).
!!
!! Written by Elizabeth Gregorio (EG)
!!

subroutine mph_iblset(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder,ibd)

#include "Flash.h"
#include "constants.h"
#include "SolidMechanics.h"

  ! Modules Used
  use SolidMechanics_data, only: sm_bodyInfo
  use sm_element_interface, only: sm_el02_mapParticles, sm_el10_mapParticles, &
                                  sm_el01_mapParticles
  use gr_sbData, ONLY : gr_sbBodyInfo
  use Driver_interface, ONLY : Driver_abortFlash

  use Grid_interface, ONLY : Grid_getListOfBlocks,   &
                             Grid_getDeltas,         &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_fillGuardCells,    &
                             Grid_getBlkBoundBox,Grid_getBlkCenterCoords

  use ib_interface, ONLY : ib_stencils

  use ImBound_data, only : ib_stencil

  implicit none

  ! IO variables
  integer, intent(in) :: sweepOrder
  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList
  real,    INTENT(IN) :: timeEndAdv,dt,dtOld
  integer, intent(in) :: ibd

!  integer, intent(in) :: blockId
!  integer, INTENT(INOUT) :: blockCount
!  integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList

  ! Internal Variables
  integer :: numPart, e, ptelem, max_ptelem, nel, p, max_ptm1
  real, allocatable, dimension(:) :: xpos,ypos
  integer, allocatable, dimension(:) :: loc_num

  ! Arugments List
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox

  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS), isAttached

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  integer :: lb,ii,jj,kk,ierr,i,j,k,dir,blockID,ind,mva,mvd

  real bsize(MDIM),coord(MDIM)

  real del(MDIM),xcell,ycell,zcell

  integer :: listofBlocks(MAXBLOCKS)
  integer :: count
  integer :: intval

  ! For the algorithm
  real, allocatable, dimension(:) :: PA, PB, P1, P0, P2, v1, v2
  real, allocatable, dimension(:) :: angl, dist, ones
  real :: u,dot,det
  integer :: nelm=2 ! Dimension for the points, 2 for (x,y) in 2-D

  nel = sm_bodyInfo(ibd)%ws_nel
  max_ptelem = maxval(sm_bodyInfo(ibd)%ws_ptelem(1:nel))
  max_ptm1 = max_ptelem-1

  allocate( xpos(max_ptelem) )
  allocate( ypos(max_ptelem) )
  allocate( PA(nelm), PB(nelm), P1(nelm), P2(nelm) )
  allocate( P0(nelm), v1(nelm), v2(nelm) )
  allocate( dist(max_ptm1), angl(max_ptm1), ones(max_ptm1) )
  allocate( loc_num(max_ptelem) )

  do lb = 1,max_ptm1
      ones(lb) = 1.0
  end do

  ! Loop through all the grid points
  do lb = 1,blockCount

        ! Loop through all the blocks
        blockID = blockList(lb)

        ! Get the block Id and boundaries
        call Grid_getBlkBoundBox(blockId,boundBox)

        ! Size of the block (far side - near side)
        bsize(:) = boundBox(2,:) - boundBox(1,:)

        ! Get co-ordinates of the block:
        call Grid_getBlkCenterCoords(blockId,coord)

        ! Get the delta x and y for the block:
        call Grid_getDeltas(blockID,del)
        
        ! Get Blocks internal limits indexes:
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

            ! Point to blocks center and face vars:
            call Grid_getBlkPtr(blockID,solnData,CENTER)
            call Grid_getBlkPtr(blockID,facexData,FACEX)
            call Grid_getBlkPtr(blockID,faceyData,FACEY)
            call Grid_getBlkPtr(blockID,facezData,FACEZ)
        

        k = 1
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
         do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           
           ! x and y coordinates for the current grid cell
           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS)

           ycell = coord(JAXIS) - bsize(JAXIS)/2.0 +   &
                   real(j - NGUARD - 1)*del(JAXIS) +  &
                   0.5*del(JAXIS)

           zcell  = 0.0 

           do k=1,max_ptelem-1

             ! End points for the line segment of the IB
             ! PA is on the left and PB is on the right
               PA = (/xpos(k), ypos(k)/)
               PB = (/xpos(k+1), ypos(k+1)/)

             ! Grid cell point
               P1 = (/xcell, ycell/)

             ! Drop a normal from P1 to the line made by connecting PA PB (not the
             ! line segment)
               u = ((P1(i)-PA(i))*(PB(i)-PA(i)) + (P1(j)-PA(j))*(PB(j)-PA(j))) / &
                   ((PB(i)-PA(i))*(PB(i)-PA(i)) + (PB(j)-PA(j))*(PB(j)-PA(j)))

             ! Re-assign u if the normal hits the line segment to the left of PA or
             ! the right of PB
               if (u .lt. 0) then
                  u = 0
               else if (u .gt. 1) then
                  u = 1
               end if 

             ! Find the point on the line segment with the shortest distance to P1
             ! (If the normal hits the line outside the line segment it is
             !  reassigned to hit the closer endpoint.)
               P0 = PA + (PB - PA)*u

             ! Determine the quadrent and angle for the "normal"
             ! (If to the left or right of the line segment the vector with the 
             !  shortest distance to the line segment will not be perpendicular)
             
               if (P0(i) .eq. PA(i) .and. P0(j) .eq. P2(j)) then
                  v1 = P1 - P0
                  v2 = P0 - PB
               else
                  v1 = P1 - P0
                  v2 = P0 - PA
               end if

               dot =   v1(1)*v2(1) + v1(2)*v2(2)
               det = -(v1(1)*v2(1) - v1(2)*v2(2))
  
            
               angl(k) = atan2(det, dot)
               dist(k) = sqrt(v1(1)**2 + v2(2)**2)

             ! Loops through all points on the grid

           end do

           ind = minloc(dist,1) ! Finds the index of the minimum distance
           mvd = dist(ind)      ! Gets the minimum distance value
           mva = angl(ind)      ! Gets teh angle for the minimum distance

           if (mva .eq. 0.0) then
               solnData(LMDA_VAR,i,j,k) = 0.0
           else
               solnData(LMDA_VAR,i,j,k) = mvd*sign(1,mva)
           end if

         end do
        end do

        k = 1
        do j=2,blkLimitsGC(HIGH,JAXIS)-1
            do i=2,blkLimitsGC(HIGH,IAXIS)-1

               solnData(NMLX_VAR,i,j,k) = -((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))/&
                                      sqrt(((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))**2)

               solnData(NMLY_VAR,i,j,k) = -((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(IAXIS))/&
                                      sqrt(((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))**2)


            end do
        end do

        k = 1
        do j=2,blkLimitsGC(HIGH,JAXIS)-1
            do i=2,blkLimitsGC(HIGH,IAXIS)-1

               solnData(TNGY_VAR,i,j,k) = ((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))/&
                                     sqrt(((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))**2+&
                                          ((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))**2)

               solnData(TNGX_VAR,i,j,k) = -((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(IAXIS))/&
                                      sqrt(((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))**2)


            end do
        end do        


        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  end do

end subroutine mph_iblset
