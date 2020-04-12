!===============================================================================
!!
!! subroutine ib_dynamic_grid_normal_vector 
!!
!! extrapolating a scalar field from inside to the outside of a region defined by sd
!!
!! THIS NEEDS TO BE BROKEN DOWN INTO MULTIPLE SUBROUTINES IN ACCORDANCE
!  WITH imBound.F90
!===============================================================================
#include "constants.h"
#include "Flash.h"

      subroutine ib_dynamic_grid_normal_vector(sd,adfx,adfy,&
                                           ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
        implicit none
        !include 'mpif.h' 

        real, dimension(:,:,:), intent(inout) :: adfx,adfy
        real, dimension(:,:,:), intent(in)    :: sd
        real :: adf

        real, intent(in)    :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        !end interface header

        integer :: i,j,k

        real :: sxl, sxr, syl, syr

        !this obtains normal components (nx,ny) of interface
        !--------normal components---------------------------------
        k = 1
        do j = jy1-NGUARD+1,jy2+NGUARD-1
          do i = ix1-NGUARD+1,ix2+NGUARD-1

              sxl = sd(i-1,j,1)
              sxr = sd(i+1,j,1)
              syl = sd(i,j-1,1)
              syr = sd(i,j+1,1)
              
              adf = sqrt( ((sxr-sxl)/2./dx)**2 + ((syr-syl)/2./dy)**2 )

              adfx(i,j,k) = (sxr-sxl)/2./dx / adf
              adfy(i,j,k) = (syr-syl)/2./dy / adf

            end do
         end do
         !--------normal components---------------------------------

        end subroutine ib_dynamic_grid_normal_vector
