!=========================================================================
!=========================================================================
!=========================================================================
        subroutine mph_getSmearedProperties2D(s,dx,dy,rho1,rho2, &
                        ix1,ix2,jy1,jy2,xnorm,ynorm, &
                        smhv,smrh)   

   
        use Multiphase_data, ONLY : mph_meshMe

        implicit none

#include "Flash.h"

        integer, intent(in) :: ix1,ix2,jy1,jy2
        real, intent(in) :: dx, dy, rho1, rho2

        real, dimension(:,:,:), intent(in):: s,xnorm,ynorm
        !integer :: icrv(NXB+2*NGUARD,NYB+2*NGUARD,1)

        !real :: th,aa,xijl,xijr, &
        !        a1,a2,cri,xid,xij,xidl,xidr,yid,yidr,yidl,yij,yijl,yijr
        integer :: i,j,k
        real, parameter :: eps = 1E-13

        !integer :: iSmear
        real :: sp
        real, dimension(:,:,:), intent(inout) :: smhv,smrh
        INTEGER, parameter :: kz1 = 1

        !real, dimension(:,:,:), intent(in out) :: xnorm,ynorm,smhv,smrhinv_facex,smrhinv_facey

!        real, dimension(NXB+2*NGUARD+1,NYB+2*NGUARD,1) :: hxl
!        real, dimension(NXB+2*NGUARD,NYB+2*NGUARD+1,1) :: hyl


         sp=1.d0*min(dx,dy)
         smhv(ix1:ix2,jy1:jy2,kz1) = &
                 (1.d0 + derf(s(ix1:ix2,jy1:jy2,kz1)/sp))/2.d0

         ! changed sides

         smrh(ix1:ix2,jy1:jy2,kz1) = &
                 (rho2/rho2) + (rho1/rho2 - rho2/rho2) * smhv(ix1:ix2,jy1:jy2,kz1)


         ! changed sides
         !smrh(ix1:ix2,jy1:jy2,kz1) = &
         !        (rho1) + (rho2 - rho1) * smhv(ix1:ix2,jy1:jy2,kz1)

!         smrminv_facex(ix1:ix2,jy1:jy2,kz1) = &




!        do j = 1,ny-1
!        do i = 1,nx
!        hxl(i,j) = (1.d0 + derf((phi(i,j)+phi(i-1,j))/2./sp))/2.d0
!        end do
!        end do
!
!        do j = 1,ny
!        do i = 1,nx-1
!        hyl(i,j) = (1.d0 + derf((phi(i,j)+phi(i,j-1))/2./sp))/2.d0
!        end do
!        end do
!
!
!        do j = 1,ny-1
!        do i = 1,nx
!        !rho_hxl(i,j) = 1.d0/(rho_g + (rho_l - rho_g) * hxl(i,j))
!        !rho_hxl(i,j) = 1.d0/rho_g + (1.d0/rho_l - 1.d0/rho_g) *
!        !hxl(i,j)
!        rho_hxl(i,j) = 1.d0/rho_g + (1.d0/rho_l - 1.d0/rho_g) *hxl_l(i,j)
!        !viscosity_hxl(i,j) = viscosity_g + (viscosity_l - viscosity_g)
!        !* hxl(i,j)
!        !viscosity_hxl(i,j) = 1.d0/(1.d0/viscosity_g + (1.d0/viscosity_l- 1.d0/viscosity_g) * hxl(i,j))
!        end do
!        end do
!
!        do j = 1,ny
!        do i = 1,nx-1
!        !rho_hyl(i,j) = 1.d0/(rho_g + (rho_l - rho_g) * hyl(i,j))
!        !rho_hyl(i,j) = 1.d0/rho_g + (1.d0/rho_l - 1.d0/rho_g) *
!        !hyl(i,j)
!        rho_hyl(i,j) = 1.d0/rho_g + (1.d0/rho_l - 1.d0/rho_g) *hyl_l(i,j)
!        !viscosity_hyl(i,j) = viscosity_g + (viscosity_l - viscosity_g)
!        !* hyl(i,j)
!        !viscosity_hyl(i,j) = 1.d0/(1.d0/viscosity_g + (1.d0/viscosity_l- 1.d0/viscosity_g) * hyl(i,j))
!        end do
!        end do
!




      end subroutine mph_getSmearedProperties2D


