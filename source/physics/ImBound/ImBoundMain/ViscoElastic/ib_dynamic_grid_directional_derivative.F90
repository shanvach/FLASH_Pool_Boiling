!===============================================================================
!!
!! subroutine ib_dynamic_grid_directional_derivative 
!!
!! defines diectional derivative sn
!!
!===============================================================================
#include "Flash.h"
#include "constants.h"

      subroutine ib_dynamic_grid_directional_derivative(sd,stest,adfx,adfy,sn,&
                                           ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
        implicit none
        !include 'mpif.h' 

        real, dimension(:,:,:), intent(inout) :: sn,adfx,adfy
        real, dimension(:,:,:), intent(in)    :: sd,stest

        real, intent(in)    :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        !end interface header

        integer :: step,i,j,k

        real :: sxl, sxr, syl, syr
        real :: h_preserve
        integer, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: pfl

        h_preserve = 0.d0

        pfl = (1 - int(sign(1.0,sd)))/2
           
        k = 1
        !------define sn(diectional derivative of level set s in normal direction)------       
        do j = jy1-NGUARD+1,jy2+NGUARD-1
           do i = ix1-NGUARD+1,ix2+NGUARD-1
 
           sn(i,j,k) = (adfx(i,j,1)*1.d0/2.d0/dx*(stest(i+1,j,1)-stest(i-1,j,1)) + & 
                        adfy(i,j,1)*1.d0/2.d0/dy*(stest(i,j+1,1)-stest(i,j-1,1)))

           end do
        end do
        !------define sn(diectional derivative of level set s in normal direction)------ 

      end subroutine ib_dynamic_grid_directional_derivative
