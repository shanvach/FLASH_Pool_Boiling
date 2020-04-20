!===============================================================================
!!
!! subroutine ib_dynamic_grid_retain_inside 
!!
!===============================================================================



      subroutine ib_dynamic_grid_retain_inside(sd,sn0,sn,&
                                           ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
        implicit none
        !include 'mpif.h' 

        real, dimension(:,:,:), intent(inout) :: sn0,sn
        real, dimension(:,:,:), intent(in)    :: sd

        real, intent(in)    :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        !end interface header

        integer :: i,j,k
        !end buit-in header

        ! retain level set inside solid. Update level set outside solid 
         do j = jy1,jy2
           do i = ix1,ix2
        if(sd(i,j,1).lt.(0.d0-dx)) then 
        sn0(i,j,1) = sn(i,j,1)
        end if
            end do
         end do
        ! retain level set inside solid. Update level set outside solid

        end subroutine ib_dynamic_grid_retain_inside
