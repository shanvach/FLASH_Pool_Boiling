!! source/physics/Multiphase/MultiphaseMain/INS/ImBound
!!
!! NAME
!!
!! mph_iblset(blockCound,blockList)
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

subroutine mph_iblset(blockCount,blockList)

#include "Flash.h"

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
  integer, intent(in) :: ibd

  ! Internal Variables
  integer :: numPart, e, ptelem, max_ptelem, nel, p
  real, allocatable, dimension(:) :: xpos,ypos,dist,angl
  real, allocatable, dimension(:) :: xacc,yacc,zacc,xnrm,ynrm,znrm
  integer, allocatable, dimension(:) :: loc_num

  ! Arugments List
  integer, intent(in) :: sweepOrder
  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox

  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS), isAttached

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  integer :: lb,blockID,ii,jj,kk,ierr,i,j,k,dir

  real bsize(MDIM),coord(MDIM)

  real del(MDIM),xcell,ycell

  integer :: listofBlocks(MAXBLOCKS)
  integer :: count
  integer :: intval

  ! Maximum number of particles per element (in the immersed boundary):
  nel = sm_bodyInfo(ibd)%ws_nel
  max_ptelem = maxval(sm_bodyInfo(ibd)%ws_ptelem(1:nel))

  allocate( xpos(max_ptelem), xnrm(max_ptelem) )
  allocate( ypos(max_ptelem), ynrm(max_ptelem) )
  allocate( dist(max_ptelem-1), angl(max_ptelm-1) )
  allocate( loc_num(max_ptelem) )

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

           ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS)  +  &
                   0.5*del(JAXIS) 

           do k=1,max_ptelem-1

           ! The following commented commands are written in a hybrid of Matlab
           ! and Fortran, they need to be ALL in Fortran and for THIS Flash code
           ! before they are uncommented.

           ! End points for the line segment of the IB
           ! PA is on the left and PB is on the right
           ! PA = [xpos(k), ypos(k)]
           ! PB = [xpos(k+1), ypos(k+1)]

           ! Grid cell point
           ! P1 = [xcell, ycell]

           ! Drop a normal from P1 to the line of PA PB
           ! u = ((P1(1)-PA(1))*(PB(1)-PA(1)) + (P1(2)-PA(2))*(PB(2)-PA(2))) / &
           !     ((PB(1)-PA(1))*(PB(1)-PA(1)) + (PB(2)-PA(2))*(PB(2)-PA(2)))

           ! Re-assign u if P1 is left of PA or right of PB
           ! if (u .lt. 0) then
           !    u = 0
           ! else if (u .gt. 1) then
           !    u = 1
           ! end if 

           ! Find the point on the line segment with the shortest distance to P1
           ! P0 = PA + (PB - PA)*u

           ! Determine the quadrent and angle for the "normal"
           ! (If to the left or right of the line segment the vector with the 
           !  shortest distance to the line segment will not be perpendicular)
           ! 
           ! if (P0(1) .eq. PA(1) .and. P0(2) .eq. P2(2)) then
           !    v1 = P1 - P0 ! REMEMBER v1 and v2 need to be declared at the top
           !                   of this file!
           !    v2 = P0 - PB
           ! else
           !    v1 = P1 - P0
           !    v2 = P0 - PA
           ! end if

           ! dot =   v1(1)*v2(1) + v1(2)*v2(2)
           ! det = -(v1(1)*v2(1) - v1(2)*v2(2))
           ! angle(k) = atan2(det, dot)

           ! d(k) = sqrt(v1(1)^2 + v2(2)^2)

           ! Loops through all points on the grid

           end do
           ! ind = find(d==min(d)) ! Would need to declare ind in the beginning!
           ! lambda(i,j) = d(ind(1))*sign(angle(ind(1)))
         end do
        end do

        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  end do

end subroutine mph_iblset
