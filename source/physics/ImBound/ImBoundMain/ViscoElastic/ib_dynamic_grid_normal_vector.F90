!===============================================================================
!!
!! subroutine ib_dynamic_grid_normal_vector 
!!
!! extrapolating a scalar field from inside to the outside of a region defined by sd
!!
!! THIS NEEDS TO BE BROKEN DOWN INTO MULTIPLE SUBROUTINES IN ACCORDANCE
!  WITH imBound.F90
!===============================================================================



      subroutine ib_dynamic_grid_normal_vector(sd,stest,adf,adfx,adfy,&
                                           ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
        implicit none
        !include 'mpif.h' 

        real, dimension(:,:,:), intent(inout) :: adfx,adfy,adf
        real, dimension(:,:,:), intent(in)    :: sd,stest

        real, intent(in)    :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        !end interface header

        integer :: i,j,k

        real :: sxl, sxr, syl, syr
        !end buit-in header

        !**args
        adfx = 0.d0
        adfy = 0.d0

        !created in this function, now passed as args
        adf = 0.d0

        !this obtains normal components (nx,ny) of interface
        !--------normal components---------------------------------
        k = 1
          !do j = 2, ny-1
          do j = jy1,jy2
            !do i = 2, nx-1
            do i = ix1,ix2

              sxl = sd(i-1,j,1)
              sxr = sd(i+1,j,1)
              syl = sd(i,j-1,1)
              syr = sd(i,j+1,1)
              
              adf(i,j,k) = sqrt( ((sxr-sxl)/2./dx)**2 + ((syr-syl)/2./dy)**2 )

              adfx(i,j,k) = (sxr-sxl)/2./dx / adf(i,j,k)
              adfy(i,j,k) = (syr-syl)/2./dy / adf(i,j,k)

            end do
         end do
         !--------normal components---------------------------------


        end subroutine ib_dynamic_grid_normal_vector
