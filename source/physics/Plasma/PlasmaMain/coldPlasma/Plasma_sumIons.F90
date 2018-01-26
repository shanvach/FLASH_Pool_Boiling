subroutine Plasma_sumIons(N_is, N_it, ix1, ix2, jy1, jy2 )

   implicit none

   real, dimension(:,:,:), intent(in) :: N_is
   real, dimension(:,:,:), intent(inout) :: N_it

   integer, intent(in) :: ix1, ix2, jy1, jy2
   integer :: i,j

   do j=jy1,jy2
     do i=ix1,ix2
         N_it(i,j,1) = N_it(i,j,1) + N_is(i,j,1)
     end do
   end do

end subroutine Plasma_sumIons
