subroutine Plasma_netCharge(dqn, dnit_pos, dnit_neg, dele, ix1, ix2, jy1, jy2) 

   !use as source term for poisson solver
   use Plasma_data, only: pls_Ce, pls_epsilon0

   implicit none

   real, dimension(:,:,:), intent(in) :: dnit_pos, dnit_neg, dele
   real, dimension(:,:,:), intent(inout) :: dqn
   
   integer, intent(in) :: ix1, ix2, jy1, jy2
   integer :: i,j

  do j=jy1,jy2
        do i=ix1,ix2
           dqn(i,j,1) = (pls_Ce/pls_epsilon0)*(dnit_pos(i,j,1) - dele(i,j,1))
        end do
  end do
 
end subroutine Plasma_netCharge
