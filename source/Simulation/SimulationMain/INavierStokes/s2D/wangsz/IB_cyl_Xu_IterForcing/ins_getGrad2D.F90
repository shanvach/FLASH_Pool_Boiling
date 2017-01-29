
  !***************************************************************
! Shizhao Wang test the pressure gradient for viscous force
! Oct 07, 2014
! Called in the ins_ab2rk3
!    call ins_getGrad2d(facexData(VELC_FACE_VAR,:,:,:),  &
!                       faceyData(VELC_FACE_VAR,:,:,:),  &
!                       slonData(PRES_VAR,:,:,:),        &
!            blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
!            blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
!             del(DIR_X),del(DIR_Y),pppxi, pppyi, pupxi, &
!                         pupyi, pvpxi, pvpyi)
  !***************************************************************

#include "constants.h"
#include "Flash.h"

      SUBROUTINE ins_getGrad2D(blockID,uni,vni,pni,ix1,ix2,jy1,jy2,dx,dy, &
           &                   ppx,ppy,pux,puy,pvx,pvy)

  !***************************************************************
  ! This subroutine computes the gradients of pressure and velocity 
  !
  ! Input:  uni,vni     = velocity at timestep n
  !         pni         = pressure at timestep n
  !         ix1,ix2     = starting and ending x indices
  !         jy1,jy2     = starting and ending y indices
  !         dx,dy       = grid spacing in x and y directions
  !
  ! Output: ppx,ppy    = gradient of pressure in x and y directions
  !         pux,puy    = gradient of u in x and y directions
  !         pvx,pvy    = gradient of v in x and y directions
  !**************************************************************

       use Grid_interface, ONLY : &
         Grid_getBlkBoundBox, Grid_getBlkCenterCoords 

      implicit none
      INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2, blockID
      REAL, INTENT(IN):: dx, dy
      REAL, DIMENSION(:,:,:), INTENT(IN):: uni, vni, pni
      REAL, DIMENSION(:,:,:), INTENT(OUT):: ppx, ppy
      REAL, DIMENSION(:,:,:), INTENT(OUT):: pux, puy, pvx, pvy

      INTEGER:: i, j
      REAL:: dx1, dy1
      ! x-component variables
      INTEGER, parameter :: kz1 = 1

      ! debug
!      real, dimension(MDIM)  :: coord,bsize
!      real ::  boundBox(2,MDIM)
!      real xedge(NXB+5),xcell(NXB+5)
!      real yedge(NYB+5),ycell(NYB+5)
!      real intsx(NXB+5), intsy(NYB+5)

      ! debug
!        intsx    = (/ (real(i), i=-2,NXB+2) /)
!        intsy    = (/ (real(i), i=-2,NYB+2) /)
        ! Get Coord and Bsize for the block:
        ! Bounding box:
!        call Grid_getBlkBoundBox(blockId,boundBox)
!        bsize(:) = boundBox(2,:) - boundBox(1,:)
        
!        call Grid_getBlkCenterCoords(blockId,coord)

        ! Point to blocks center and face vars:
!        xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
!        xcell = xedge(:) + dx/2.0;
 
!        yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
!        ycell = yedge(:) + dy/2.0;


!        print*, 'boundBox',  boundBox
!        print*, 'bsize', bsize
!        print*, 'coord', coord
!        print*, 'xcell', xcell(1), xcell(NXB+1)
!        print*, 'ycell', ycell(1), ycell(NYB+1)


      ! grid spacings
      dx1 = 0.5/dx
      dy1 = 0.5/dy

      !++++++++++  Pressure   ++++++++++
       do j = jy1,jy2
          do i = ix1,ix2
             ppx(i,j,kz1) = (pni(i+1,j,kz1) - pni(i-1,j,kz1))*dx1 
             ppy(i,j,kz1) = (pni(i,j+1,kz1) - pni(i,j-1,kz1))*dy1 
!            ppx(i,j,kz1) = xcell(i) + ycell(j)*ycell(j)
!             ppy(i,j,kz1) = sin(xcell(i)) + &
!                            cos(ycell(j))  
          enddo
       enddo
      
      !++++++++++  U-COMPONENT  ++++++++++
       do j = jy1,jy2
          do i = ix1,ix2
             pux(i,j,kz1) = (uni(i+1,j,kz1) - uni(i-1,j,kz1))*dx1 
             puy(i,j,kz1) = (uni(i,j+1,kz1) - uni(i,j-1,kz1))*dy1 
          enddo
       enddo

      !++++++++++  V-COMPONENT  ++++++++++
       do j = jy1,jy2
          do i = ix1,ix2
             pvx(i,j,kz1) = (vni(i+1,j,kz1) - vni(i-1,j,kz1))*dx1 
             pvy(i,j,kz1) = (vni(i,j+1,kz1) - vni(i,j-1,kz1))*dy1 
          enddo
       enddo


       RETURN
       END SUBROUTINE ins_getGrad2D


