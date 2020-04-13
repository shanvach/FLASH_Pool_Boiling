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

        real :: sxplus, sxmins, syplus, symins, up, vp
        real :: h_preserve
        integer, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: pfl

        h_preserve = 0.d0

        pfl = (1 - int(sign(1.0,sd)))/2
      
        sn = 0.0d0
     
        k = 1
        !------define sn(diectional derivative of level set s in normal direction)------       
        do j = jy1,jy2
           do i = ix1,ix2
 
          !normal vectors
          up = adfx(i,j,k)
          vp = adfy(i,j,k)

          !gradients
          sxplus = (stest(i+1,j,k) - stest(i,j,k)) / dx
          sxmins = (stest(i,j,k) - stest(i-1,j,k)) / dx
          syplus = (stest(i,j+1,k) - stest(i,j,k)) / dy
          symins = (stest(i,j,k) - stest(i,j-1,k)) / dy

          ! use dx/2 as dt to advect level set
          sn(i,j,k) = (max(up,0.0d0)*sxmins + min(up,0.0d0)*sxplus + &
                       max(vp,0.0d0)*symins + min(vp,0.0d0)*syplus)*pfl(i,j,k)

           end do
        end do
        !------define sn(diectional derivative of level set s in normal direction)------ 

      end subroutine ib_dynamic_grid_directional_derivative
