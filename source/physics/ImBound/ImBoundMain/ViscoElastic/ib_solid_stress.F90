!===============================================================================
!!
!! subroutine ib_solid_stress
!!
!! calculate deformation gradient tensor , solid strain tensor and strain magnitude 
!!
!! args inout:: Tau       - strain tensor (2D=4elements)
!! args in   :: sd        - level set for phase interface
!!           :: sX,sY     - dynamic grid (X,Y)
!!           :: A,AT,A_inv- deformation gradient tensor, transpose, and
!!                          inverse
!!           :: Taux,Tauy - Div(Tau),divergence of tensor 
!!           :: Taum      - magnitude(Taux,Tauy)
!!           :: dx,dy     - spacing in x and y directions
!!           :: ix1,ix2   - low,high block x indeces
!!           :: jy1,jy2   - low,high block y indeces 
!===============================================================================
#include "constants.h"
#include "Flash.h"

        subroutine ib_solid_stress(sd,sX,sY,Tau1,Tau2,Tau3,Tau4,&
                                   ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)   
        implicit none
        !include 'mpif.h'
        real, dimension(:,:,:), intent(inout) :: Tau1, Tau2, Tau3, Tau4
        real, dimension(:,:,:), intent(inout)    :: sd,sX,sY
        real :: A1,A2,A3,A4
        real :: AT1,AT2,AT3,AT4
        real :: A_inv1,A_inv2,A_inv3,A_inv4
        !real, dimension(:,:,:)                :: Taux,Tauy,Taum

        real, intent(in)    :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2

        integer :: i,j,k
        real    :: ul,ur,vl,vr

        Tau1 = 0.d0
        Tau2 = 0.d0
        Tau3 = 0.d0
        Tau4 = 0.d0

        !Taux = 0.d0
        !Tauy = 0.d0
        !Taum = 0.d0 

        !A(i,j,ind) is the deformation gradient tensor in 2D 
        !ind is the index of each component in the deformation gradient tensor. dX/dx:1; dX/dy:2; dY/dx:3; dY/dy:4.
        k = 1
        do j = jy1-NGUARD+1,jy2+NGUARD-1
          do i = ix1-NGUARD+1,ix2+NGUARD-1

          A1     = 0.d0
          A2     = 0.d0
          A3     = 0.d0
          A4     = 0.d0

          AT1    = 0.d0
          AT2    = 0.d0
          AT3    = 0.d0
          AT4    = 0.d0

          A_inv1 = 1.d0/(2*dx)*(sX(i+1,j,k)-sX(i-1,j,k))
          A_inv2 = 1.d0/(2*dy)*(sX(i,j+1,k)-sX(i,j-1,k))
          A_inv3 = 1.d0/(2*dx)*(sY(i+1,j,k)-sY(i-1,j,k))
          A_inv4 = 1.d0/(2*dy)*(sY(i,j+1,k)-sY(i,j-1,k))


          if(abs(A_inv1*A_inv4-A_inv2*A_inv3).gt.1E-12) then

           A1 = 1.d0/(A_inv1*A_inv4-A_inv2*A_inv3)*A_inv4

           A2 = 1.d0/(A_inv1*A_inv4-A_inv2*A_inv3)*(-A_inv2) 

           A3 = 1.d0/(A_inv1*A_inv4-A_inv2*A_inv3)*(-A_inv3) 

           A4 = 1.d0/(A_inv1*A_inv4-A_inv2*A_inv3)*A_inv1

          end if

          AT1 = A1
          AT2 = A3
          AT3 = A2
          AT4 = A4

          !Tau is the strain tensor
          Tau1(i,j,k) = (A1*AT1+A2*AT3-1)
          Tau2(i,j,k) = (A1*AT2+A2*AT4)
          Tau3(i,j,k) = (A3*AT1+A4*AT3)
          Tau4(i,j,k) = (A3*AT2+A4*AT4-1)
        
          end do
        end do

       ! !divergence of Tau gives solid stress in x and y directions
       ! k = 1
       ! do j = jy1,jy2
       !   do i = ix1,ix2
       !   !if(sd(i,j,1).le.0.d0) then
       !   Taux(i,j,k) = 1.d0/(2*dx)*(Tau(i+1,j,1)-Tau(i-1,j,1)) + 1.d0/(2*dy)*(Tau(i,j+1,2)-Tau(i,j-1,2)) 
       !   Tauy(i,j,k) = 1.d0/(2*dx)*(Tau(i+1,j,3)-Tau(i-1,j,3)) + 1.d0/(2*dy)*(Tau(i,j+1,4)-Tau(i,j-1,4))
       !   Taum(i,j,k) = sqrt(Taux(i,j,k)**2 + Tauy(i,j,k)**2)
       ! !end if 
       !    end do
       ! end do

      end subroutine ib_solid_stress
