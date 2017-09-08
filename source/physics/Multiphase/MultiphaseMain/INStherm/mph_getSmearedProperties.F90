subroutine mph_getSmearedProperties2D(s,pf,dx,dy,rho1,rho2,ix1,ix2,jy1,jy2,xnorm,ynorm,smhv,smrh)   

   
        use Multiphase_data, ONLY : mph_meshMe

        implicit none

#include "Flash.h"

        integer, intent(in) :: ix1,ix2,jy1,jy2
        real, intent(in) :: dx, dy, rho1, rho2

        real, dimension(:,:,:), intent(in):: s,pf,xnorm,ynorm

        integer :: i,j,k
        real, parameter :: eps = 1E-13

        real :: sp,pi
        real, dimension(:,:,:), intent(inout) :: smhv,smrh
        INTEGER, parameter :: kz1 = 1

         pi = acos(-1.0)

         !sp=0.5*dx
         !smhv(ix1:ix2,jy1:jy2,kz1) = &
         !        (1.d0 + derf(s(ix1:ix2,jy1:jy2,kz1)/sp))/2.d0

         !smrh(ix1:ix2,jy1:jy2,kz1) = rho2/rho2 + (rho2/rho1 - rho2/rho2)*smhv(ix1:ix2,jy1:jy2,kz1)

         sp = 1.5*dx

         do j=jy1,jy2
           do i=ix1,ix2

              !if(abs(s(i,j,kz1)) .le. sp) then ! Symmetric smearing - AD
              !!if(abs(s(i,j,kz1)) .le. sp .and. s(i,j,kz1) .lt. 0.0) then ! Asymmetric smearing - AD

              !smhv(i,j,kz1) = 0.5 + s(i,j,kz1)/(2*sp) + sin(2*pi*s(i,j,kz1)/(2*sp))/(2*pi)

              !else

                  if(s(i,j,kz1) .ge. 0.0) then

                        smhv(i,j,kz1) = 1.0
  
                  else

                        smhv(i,j,kz1) = 0.0

                  end if

              !end if


           end do
         end do

         smrh(ix1:ix2,jy1:jy2,kz1) = &
                 (rho2/rho2) + (rho1/rho2 - rho2/rho2) * smhv(ix1:ix2,jy1:jy2,kz1)

         smrh(ix1:ix2,jy1:jy2,kz1) = 1./smrh(ix1:ix2,jy1:jy2,kz1)

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

        real :: sp,pi
        real, dimension(:,:,:), intent(inout) :: smhv,smrh


       pi = acos(-1.0)

       !sp=0.5*dx
       !smhv(ix1:ix2,jy1:jy2,kz1:kz2) = &
       !         (1.d0 + derf(s(ix1:ix2,jy1:jy2,kz1:kz2)/sp))/2.d0

       !smrh(ix1:ix2,jy1:jy2,kz1:kz2) = rho2/rho2 + (rho2/rho1 - rho2/rho2)*smhv(ix1:ix2,jy1:jy2,kz1:kz2)

       sp = 1.5*dx

       do k=kz1,kz2
        do j=jy1,jy2
          do i=ix1,ix2

              !if(abs(s(i,j,k)) .le. sp) then ! Symmetric smearing - AD
              !!if(abs(s(i,j,k)) .le. sp .and. s(i,j,k) .lt. 0.0) then ! Asymmetric smearing - AD

              !smhv(i,j,k) = 0.5 + s(i,j,k)/(2*sp) + sin(2*pi*s(i,j,k)/(2*sp))/(2*pi)

              !else

                  if(s(i,j,k) .ge. 0.0) then

                        smhv(i,j,k) = 1.0
  
                  else

                        smhv(i,j,k) = 0.0

                  end if

              !end if


           end do
        end do
       end do      

       smrh(ix1:ix2,jy1:jy2,kz1:kz2) = &
               (rho2/rho2) + (rho1/rho2 - rho2/rho2) * smhv(ix1:ix2,jy1:jy2,kz1:kz2)

       smrh(ix1:ix2,jy1:jy2,kz1:kz2) = 1./smrh(ix1:ix2,jy1:jy2,kz1:kz2)


end subroutine mph_getSmearedProperties3D


