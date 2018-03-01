subroutine Heat_Solve(T_p, T_o, u, v, dt, dx, dy, dz, inRe, ix1,ix2, jy1, jy2, T_res)

  use Heat_AD_data

#include "Heat_AD.h"

  implicit none
  real, dimension(:,:,:), intent(inout) :: T_p
  real, dimension(:,:,:), intent(in) :: T_o
  real, dimension(:,:,:), intent(in) :: u,v
  real, intent(in) :: dt, dx, dy, dz, inRe
  integer, intent(in) :: ix1, ix2, jy1, jy2

  real, intent(out) :: T_res

  integer :: i,j

  real :: u_plus, u_mins, v_plus, v_mins, u_conv, v_conv
  real :: Tx_plus, Tx_mins, Ty_plus, Ty_mins

#ifdef HEAT_CENTRAL

  ! with Loop Vectorization ! (Second Order Central)

  T_p(ix1:ix2,jy1:jy2,1) = T_o(ix1:ix2,jy1:jy2,1) &
  +((dt*inRe)/(ht_Pr*dx*dx))*(T_o(ix1+1:ix2+1,jy1:jy2,1)+T_o(ix1-1:ix2-1,jy1:jy2,1)-2*T_o(ix1:ix2,jy1:jy2,1))&
  +((dt*inRe)/(ht_Pr*dy*dy))*(T_o(ix1:ix2,jy1+1:jy2+1,1)+T_o(ix1:ix2,jy1-1:jy2-1,1)-2*T_o(ix1:ix2,jy1:jy2,1))&
  -((dt*(u(ix1+1:ix2+1,jy1:jy2,1) + u(ix1:ix2,jy1:jy2,1))/2)/(2*dx))*(T_o(ix1+1:ix2+1,jy1:jy2,1)-T_o(ix1-1:ix2-1,jy1:jy2,1))&
  -((dt*(v(ix1:ix2,jy1+1:jy2+1,1) + v(ix1:ix2,jy1:jy2,1))/2)/(2*dy))*(T_o(ix1:ix2,jy1+1:jy2+1,1)-T_o(ix1:ix2,jy1-1:jy2-1,1))

#endif


#ifdef HEAT_UPWIND
  ! without Loop Vectorization ! (First Order Upwind)

  do j=jy1,jy2
     do i=ix1,ix2

     u_conv = (u(i+1,j,1)+u(i,j,1))/2.
     v_conv = (v(i,j+1,1)+v(i,j,1))/2.

     u_plus = max(u_conv, 0.)
     u_mins = min(u_conv, 0.)

     v_plus = max(v_conv, 0.)
     v_mins = min(v_conv, 0.)

     Tx_plus = T_o(i+1,j,1)-T_o(i,j,1)
     Tx_mins = T_o(i,j,1)-T_o(i-1,j,1)

     Ty_plus = T_o(i,j+1,1)-T_o(i,j,1)
     Ty_mins = T_o(i,j,1)-T_o(i,j-1,1)

     T_p(i,j,1) = T_o(i,j,1) + ((dt*inRe)/(ht_Pr*dx*dx))*(T_o(i+1,j,1)+T_o(i-1,j,1)-2.*T_o(i,j,1))&
                             + ((dt*inRe)/(ht_Pr*dy*dy))*(T_o(i,j+1,1)+T_o(i,j-1,1)-2.*T_o(i,j,1))&
                             - ((dt)/dx) * (u_plus*Tx_mins + u_mins*Tx_plus)&
                             - ((dt)/dy) * (v_plus*Ty_mins + v_mins*Ty_plus)

     end do
  end do 

#endif

  do i = ix1,ix2
     T_res = T_res + sum((T_o(i,:,1)-T_p(i,:,1))**2)
  end do

  T_res = sqrt(T_res/size(T_o))

  !print *,T_res

end subroutine Heat_Solve
