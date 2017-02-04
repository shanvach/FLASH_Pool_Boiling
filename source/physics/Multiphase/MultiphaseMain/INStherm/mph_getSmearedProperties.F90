subroutine mph_getSmearedProperties2D(s,pf,dx,dy,rho1,rho2,ix1,ix2,jy1,jy2,xnorm,ynorm,smhv,smrh)   

   
        use Multiphase_data, ONLY : mph_meshMe

        implicit none

#include "Flash.h"

        integer, intent(in) :: ix1,ix2,jy1,jy2
        real, intent(in) :: dx, dy, rho1, rho2

        real, dimension(:,:,:), intent(in):: s,pf,xnorm,ynorm

        integer :: i,j,k
        real, parameter :: eps = 1E-13

        real :: sp
        real, dimension(:,:,:), intent(inout) :: smhv,smrh
        INTEGER, parameter :: kz1 = 1

         !sp=1.d0*min(dx,dy)
         !smhv(ix1:ix2,jy1:jy2,kz1) = &
         !        (1.d0 + derf(s(ix1:ix2,jy1:jy2,kz1)/sp))/2.d0

         smhv(ix1:ix2,jy1:jy2,kz1) = pf(ix1:ix2,jy1:jy2,kz1)

         smrh(ix1:ix2,jy1:jy2,kz1) = &
                 (rho2/rho2) + (rho2/rho1 - rho2/rho2) * smhv(ix1:ix2,jy1:jy2,kz1)

end subroutine mph_getSmearedProperties2D

subroutine mph_getSmearedProperties3D(s,pf,dx,dy,dz,rho1,rho2,ix1,ix2,jy1,jy2,kz1,kz2,xnorm,ynorm,znorm,smhv,smrh) 

   
        use Multiphase_data, ONLY : mph_meshMe

        implicit none

#include "Flash.h"

        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        real, intent(in) :: dx, dy, dz, rho1, rho2

        real, dimension(:,:,:), intent(in):: s,pf,xnorm,ynorm,znorm

        integer :: i,j,k
        real, parameter :: eps = 1E-13

        real :: sp
        real, dimension(:,:,:), intent(inout) :: smhv,smrh

         !sp=1.d0*min(dx,dy)
         !smhv(ix1:ix2,jy1:jy2,kz1:kz2) = &
         !        (1.d0 + derf(s(ix1:ix2,jy1:jy2,kz1:kz2)/sp))/2.d0

         smhv(ix1:ix2,jy1:jy2,kz1:kz2) = pf(ix1:ix2,jy1:jy2,kz1:kz2)

         smrh(ix1:ix2,jy1:jy2,kz1:kz2) = &
                 (rho2/rho2) + (rho2/rho1 - rho2/rho2) * smhv(ix1:ix2,jy1:jy2,kz1:kz2)

end subroutine mph_getSmearedProperties3D


