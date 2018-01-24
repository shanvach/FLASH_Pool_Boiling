subroutine Plasma_sumNeutrals(N_as, N_at, ix1, ix2, jy1, jy2)

   implicit none

   real, dimension(:,:,:), intent(in) :: N_as
   real, dimension(:,:,:), intent(inout) :: N_at
   
   integer, intent(in) :: ix1, ix2, jy1, jy2
   integer :: i,j

   do j=jy1,jy2
     do i=ix1,ix2
         N_at(i,j,1) = N_at(i,j,1) + N_as(i,j,1)
     end do
   end do

end subroutine Plasma_sumNeutrals
