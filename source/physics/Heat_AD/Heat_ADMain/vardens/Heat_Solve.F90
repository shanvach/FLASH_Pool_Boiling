subroutine Heat_Solve(T_p,T_o,T_rhs,dt,ix1,ix2,jy1,jy2,kz1,kz2,T_res)

  use Heat_AD_data

#include "Heat_AD.h"

  implicit none
  real, dimension(:,:,:), intent(inout) :: T_p
  real, dimension(:,:,:), intent(in) :: T_o
  real, dimension(:,:,:), intent(in) :: T_rhs
  real, intent(in) :: dt
  integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
 
  real, intent(inout) :: T_res
  integer :: i,j,k

  !do j=jy1,jy2
  ! do i=ix1,ix2
  !     T_p(i,j,1) = T_o(i,j,1) + dt * T_rhs(i,j,1)
  ! end do
  !end do

  T_p(ix1:ix2,jy1:jy2,kz1:kz2) = T_o(ix1:ix2,jy1:jy2,kz1:kz2) + dt*T_rhs(ix1:ix2,jy1:jy2,kz1:kz2)

  T_res = sum(sum(sum((T_o(:,:,:)-T_p(:,:,:))**2,1),1))

  !do i = ix1,ix2
  !   T_res = T_res + sum((T_o(i,:,:)-T_p(i,:,:))**2)
  !end do

  T_res = sqrt(T_res/size(T_o))

  !print *,T_res

end subroutine Heat_Solve
