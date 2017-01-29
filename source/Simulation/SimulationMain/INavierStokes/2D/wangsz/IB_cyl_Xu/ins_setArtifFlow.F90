
  !***************************************************************
! Shizhao Wang set the artifical flows
! Oct 07, 2014
! Called in the ins_ab2rk3
!    call ins_setArtifFlow(facexData(VELC_FACE_VAR,:,:,:),  &
!                       faceyData(VELC_FACE_VAR,:,:,:),  &
!                       slonData(PRES_VAR,:,:,:))
  !***************************************************************

#include "constants.h"
#include "Flash.h"

      SUBROUTINE ins_setArtifFlow(blockID,uni,vni,pni,ix1,ix2,jy1,jy2,dx,dy)

  !***************************************************************
  ! This subroutine computes the gradients of pressure and velocity 
  !
  ! Input:  uni,vni     = velocity at timestep n
  !         pni         = pressure at timestep n
  !         ix1,ix2     = starting and ending x indices
  !         jy1,jy2     = starting and ending y indices
  !         dx,dy       = grid spacing in x and y directions
  !
  !**************************************************************

       use Grid_interface, ONLY : &
         Grid_getBlkBoundBox, Grid_getBlkCenterCoords 

      implicit none
      INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2, blockID
      REAL, INTENT(IN):: dx, dy
      REAL, DIMENSION(:,:,:), INTENT(OUT) :: uni, vni, pni

      INTEGER:: i, j
      REAL:: dx1, dy1
      ! x-component variables
      INTEGER, parameter :: kz1 = 1

      ! debug
      real, dimension(MDIM)  :: coord,bsize
      real ::  boundBox(2,MDIM)
      real xedge(NXB+5),xcell(NXB+5)
      real yedge(NYB+5),ycell(NYB+5)
      real intsx(NXB+5), intsy(NYB+5)
      real pi, x, y

      ! debug
        pi = acos(-1.)

        intsx    = (/ (real(i), i=-2,NXB+2) /)
        intsy    = (/ (real(i), i=-2,NYB+2) /)
        ! Get Coord and Bsize for the block:
        ! Bounding box:
        call Grid_getBlkBoundBox(blockId,boundBox)
        bsize(:) = boundBox(2,:) - boundBox(1,:)
        
        call Grid_getBlkCenterCoords(blockId,coord)

        ! Point to blocks center and face vars:
        xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + dx*intsx;
        xcell = xedge(:) + dx/2.0;
 
        yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + dy*intsy;
        ycell = yedge(:) + dy/2.0;


!        print*, 'boundBox',  boundBox
!        print*, 'bsize', bsize
!        print*, 'coord', coord
!        print*, 'xcell', xcell(1), xcell(NXB+1)
!        print*, 'ycell', ycell(1), ycell(NYB+1)


      ! grid spacings
      dx1 = 2.0/dx
      dy1 = 2.0/dy

      !++++++++++  Pressure   ++++++++++
       do j = jy1,jy2
          y = ycell(j)
          do i = ix1,ix2
             x = xcell(i)
             pni(i,j,kz1) = -0.25*(cos(2.*pi*x) + cos(2.*pi*y) ) 
          enddo
       enddo
      
      !++++++++++  U-COMPONENT  ++++++++++
       do j = jy1,jy2
          y = ycell(j)
          do i = ix1,ix2+1
             x = xedge(i)
             uni(i,j,kz1) = -cos(pi*x)*sin(pi*y)
          enddo
       enddo

      !++++++++++  V-COMPONENT  ++++++++++
       do j = jy1,jy2+1
          y = yedge(j)
          do i = ix1,ix2
             x = xcell(i)
             vni(i,j,kz1) = sin(pi*x)*cos(pi*y)
          enddo
       enddo


       RETURN
       END SUBROUTINE ins_setArtifFlow


