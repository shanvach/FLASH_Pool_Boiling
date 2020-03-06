!===============================================================================
!!
!! subroutine ib_dynamic_grid_directional_derivative 
!!
!! defines diectional derivative sn
!!
!===============================================================================



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
        !end buit-in header

        sn = 0.d0

        !created in this function, now passed as args
        !sn0 = 0.d0
        ! define thickness of the region outside the interface that is not extrapolated(eg. diffused region)
        ! h_preserve = 0 if only preserve level set inside the interface(sd=0)
        h_preserve = 0.d0

       
        !------define sn(diectional derivative of level set s in normal direction)------
           !do j = 2, ny-1 
           do j = jy1,jy2
              !do i = 2, nx-1
              do i = ix1,ix2
              if( sd(i,j,1).lt.(h_preserve-dx) ) then
              sn(i,j,1) = adfx(i,j,1)*1.d0/2.d0/dx*(stest(i+1,j,1)-stest(i-1,j,1)) + & 
                          adfy(i,j,1)*1.d0/2.d0/dy*(stest(i,j+1,1)-stest(i,j-1,1))
              end if
              end do
           end do
        !------define sn(diectional derivative of level set s in normal direction)------ 



      end subroutine ib_dynamic_grid_directional_derivative
