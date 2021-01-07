!!===============================================================================
!!
!! subroutine to assign density, viscoity and scalars on grid
!!
!=============================================================================== 

      subroutine ib_fluid_props(s,rho1, rho2, xmu1, xmu2, xmus,&
                                rho, xmu, xms, blockID,     &
                                ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)

#include "constants.h"
#include "Flash.h"

        use Grid_interface, only: Grid_getBlkBoundBox, Grid_getBlkCenterCoords, Grid_getDeltas
        use Grid_data, only: gr_meshMe

        implicit none

        include "Flash_mpi.h"

        real, dimension(:,:,:), intent(inout) :: s, rho, xmu, xms
        real, intent(in)    :: rho1, rho2, xmu1, xmu2, xmus
        real, intent(in)    :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        integer, intent(in) :: blockID

        real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: phi, s1, dns

        integer :: nx, ny, nz
        integer :: i,j,k,ierr
        integer :: nip,nip_max,inv,n1, n3, nn1, nn2, n, m, n2 
        real    :: ul,ur,vl,vr
        real    :: th, xx, yy, xx1, yy1, xx2, yy2
        real    :: eps, psi
        real    :: r1, r2, r3
        real    :: xx01, yy01, xx02, yy02, xx03, yy03
        real    :: xcell1,ycell1,zcell1,xcell2,ycell2,zcell2

        real, dimension((ix2-ix1+1)*(jy2-jy1+1),6)  :: xip
        real, dimension(2,MDIM) :: boundBox
        real bsize(MDIM),coord(MDIM)
        real del(MDIM)


        call Grid_getBlkBoundBox(blockID,boundBox)
        bsize(:) = boundBox(2,:) - boundBox(1,:)

        call Grid_getBlkCenterCoords(blockID,coord)

        call Grid_getDeltas(blockID,del)

        eps = 1.E-10
       
        phi = s
        do j = jy1, jy2
          do i = ix1, ix2
                 psi = (1.d0 + derf(-phi(i,j,1)/(2.d0*dx)))/2.d0   !
                 xmu(i,j,1) = psi*(xmu1-xmu2) + xmu2               !
                 !set xmu inside solid to be xmu1, outside to be xmu2
                 !set xms inside solid to be xmus, outside to be 0
                 xms(i,j,1) = psi*(xmus-0.d0) + 0.d0
                 rho(i,j,1) = psi*(rho1-rho2) + rho2               !
             end do
        end do        

      end subroutine ib_fluid_props
