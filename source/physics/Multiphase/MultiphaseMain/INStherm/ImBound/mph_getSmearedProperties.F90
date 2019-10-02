subroutine mph_getSmearedProperties2D(s,pf,dx,dy,rho1,rho2,ix1,ix2,jy1,jy2,xnorm,ynorm,smhv,smrh,lambda)   

   
        use Multiphase_data, ONLY : mph_meshMe

        implicit none

#include "Flash.h"
#include "constants.h"

        integer, intent(in) :: ix1,ix2,jy1,jy2
        real, intent(in) :: dx, dy, rho1, rho2

        real, dimension(:,:,:), intent(in):: s,pf,xnorm,ynorm,lambda

        integer :: i,j,k
        real, parameter :: eps = 1E-13

        real :: sp,pi
        real, dimension(:,:,:), intent(inout) :: smhv,smrh
        INTEGER, parameter :: kz1 = 1

        real :: sunion(NXB+2*NGUARD,NYB+2*NGUARD,1)
        real :: pfl(NXB+2*NGUARD,NYB+2*NGUARD,1)
        real :: rho3

         rho3 = (rho2 + rho1)/2.0

         pi = acos(-1.0)

         sp = 1.5*dx

         sunion = s
         pfl = 0.0

         k = kz1

         !do j=jy1,jy2
         !   do i=ix1,ix2
         !       sunion(i,j,kz1) = min(s(i,j,kz1),-lambda(i,j,kz1))
         !   end do
         !end do

         pfl(ix1-1:ix2+1,jy1-1:jy2+1,k)  = 0.0
         pfl(ix1-1:ix2+1,jy1-1:jy2+1,k)  = (sign(1.0,lambda(ix1-1:ix2+1,jy1-1:jy2+1,k))+1.0)/2.0 

         do j=jy1,jy2
           do i=ix1,ix2

              if(abs(sunion(i,j,kz1)) .le. sp .and. lambda(i,j,kz1) .lt. 0.0) then ! Symmetric smearing - AD

              smhv(i,j,kz1) = 0.5 + sunion(i,j,kz1)/(2*sp) + sin(2*pi*sunion(i,j,kz1)/(2*sp))/(2*pi)

              else

                  if(sunion(i,j,kz1) .ge. 0.0) then

                        smhv(i,j,kz1) = 1.0
  
                  else

                        smhv(i,j,kz1) = 0.0

                  end if

              end if

           end do
         end do

         do j=jy1,jy2
          do i=ix1,ix2

              smrh(i,j,k) = (1-smhv(i,j,k))*(1-pfl(i,j,k))*(rho2/rho2) + &
                              (smhv(i,j,k))*(1-pfl(i,j,k))*(rho1/rho2) + &
                             (pfl(i,j,k))*(rho3/rho2)

          end do
         end do

         smrh(ix1:ix2,jy1:jy2,kz1) = 1./smrh(ix1:ix2,jy1:jy2,kz1)

end subroutine mph_getSmearedProperties2D

subroutine mph_getSmearedProperties3D(s,pf,dx,dy,dz,rho1,rho2,ix1,ix2,jy1,jy2,kz1,kz2,xnorm,ynorm,znorm,smhv,smrh,lambda) 

   
        use Multiphase_data, ONLY : mph_meshMe

        implicit none

#include "Flash.h"
#include "constants.h"

        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        real, intent(in) :: dx, dy, dz, rho1, rho2

        real, dimension(:,:,:), intent(in):: s,pf,xnorm,ynorm,znorm,lambda

        integer :: i,j,k
        real, parameter :: eps = 1E-13

        real :: sp,pi
        real, dimension(:,:,:), intent(inout) :: smhv,smrh

        real :: sunion(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD)
        real :: pfl(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD)
        real :: rho3

       rho3 = (rho2 + rho1)/2.0

       pi = acos(-1.0)

       sp = 1.5*dx

       sunion = s
       pfl = 0.0

       !do k=kz1,kz2 
       ! do j=jy1,jy2
       !     do i=ix1,ix2
       !         sunion(i,j,k) = min(s(i,j,k),-lambda(i,j,k))
       !     end do
       !  end do
       !end do

       pfl(ix1-1:ix2+1,jy1-1:jy2+1,kz1-1:kz2+1)  = 0.0
       pfl(ix1-1:ix2+1,jy1-1:jy2+1,kz1-1:kz2+1)  = (sign(1.0,lambda(ix1-1:ix2+1,jy1-1:jy2+1,kz1-1:kz2+1))+1.0)/2.0 


       do k=kz1,kz2
        do j=jy1,jy2
          do i=ix1,ix2

              if(abs(sunion(i,j,k)) .le. sp .and. lambda(i,j,kz1) .lt. 0.0) then

              smhv(i,j,k) = 0.5 + sunion(i,j,k)/(2*sp) + sin(2*pi*sunion(i,j,k)/(2*sp))/(2*pi)

              else

                  if(sunion(i,j,k) .ge. 0.0) then

                        smhv(i,j,k) = 1.0
  
                  else

                        smhv(i,j,k) = 0.0

                  end if

              end if

           end do
        end do
       end do      

       do k=kz1,kz2
         do j=jy1,jy2
          do i=ix1,ix2

              smrh(i,j,k) = (1-smhv(i,j,k))*(1-pfl(i,j,k))*(rho2/rho2) + &
                              (smhv(i,j,k))*(1-pfl(i,j,k))*(rho1/rho2) + &
                             (pfl(i,j,k))*(rho3/rho2)

          end do
         end do
      end do

       smrh(ix1:ix2,jy1:jy2,kz1:kz2) = 1./smrh(ix1:ix2,jy1:jy2,kz1:kz2)


end subroutine mph_getSmearedProperties3D
