subroutine Plasma_Solve(T_p, T_gen, T_o, T_source, dfun, dcoeff, dt, dx, dy, ix1,ix2, jy1, jy2, T_res)

#include "Plasma.h"

  implicit none
  real, dimension(:,:,:), intent(inout) :: T_p
  real, dimension(:,:,:), intent(in) :: T_o, T_source, T_gen, dfun, dcoeff
  real, intent(in) :: dt, dx, dy
  integer, intent(in) :: ix1, ix2, jy1, jy2

  real,intent(out) :: T_res

  integer :: i,j

  do j=jy1,jy2
     do i=ix1,ix2

     !T_p(i,j,1) = T_o(i,j,1) + 0.5*(1.0-sign(1.0,dfun(i,j,1)))*((dt*dcoeff(i,j,1))/(dx*dx))*(T_o(i+1,j,1)+T_o(i-1,j,1)-2.*T_o(i,j,1))&
     !                        + 0.5*(1.0-sign(1.0,dfun(i,j,1)))*((dt*dcoeff(i,j,1))/(dy*dy))*(T_o(i,j+1,1)+T_o(i,j-1,1)-2.*T_o(i,j,1))&
     !                        + 0.5*(1.0-sign(1.0,dfun(i,j,1)))*(dt*T_gen(i,j,1))

     T_p(i,j,1) = T_o(i,j,1) + (dt*T_gen(i,j,1))   & 
                             + (dt*T_source(i,j,1)) &
                             + ((dt*dcoeff(i,j,1))/(dx*dx))*(T_o(i+1,j,1)+T_o(i-1,j,1)-2.*T_o(i,j,1)) &
                             + ((dt*dcoeff(i,j,1))/(dy*dy))*(T_o(i,j+1,1)+T_o(i,j-1,1)-2.*T_o(i,j,1))  

     !set minimum value to avoid negative density
     if (T_p(i,j,1).le.1e6) then
        T_p(i,j,1) = 1e6
     end if

     end do
  end do 

  do i = ix1,ix2
     T_res = T_res + sum((T_o(i,:,1)-T_p(i,:,1))**2)
  end do

  T_res = sqrt(T_res/size(T_o))

end subroutine Plasma_Solve
