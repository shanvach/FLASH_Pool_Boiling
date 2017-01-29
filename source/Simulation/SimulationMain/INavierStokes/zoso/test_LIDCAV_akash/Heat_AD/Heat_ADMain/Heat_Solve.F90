subroutine Heat_Solve(T_p, T_o, u, v, dt, dx, dy, dz, alfa, inRe, ix1,ix2, jy1, jy2)

  implicit none

  real, dimension(:,:,:), intent(inout) :: T_p
  real, dimension(:,:,:), intent(in) :: T_o
  real, dimension(:,:,:), intent(in) :: u,v
  real, intent(in) :: dt, dx, dy, dz, inRe
  integer, intent(in) :: ix1, ix2, jy1, jy2
  real, intent(in) :: alfa

  real :: Pr

  real :: T_res

  integer :: i,j

  real :: u_plus, u_mins, v_plus, v_mins, u_conv, v_conv
  real :: Tx_plus, Tx_mins, Ty_plus, Ty_mins
  real :: u_plus1, u_mins1, v_plus1, v_mins1
  real :: Tx, Ty

  Pr = 0.6


  ! with Loop Vectorization !

  !T_p(ix1:ix2,jy1:jy2,1) = T_o(ix1:ix2,jy1:jy2,1) &
  !+((alfa*dt*inRe)/(Pr*dx*dx))*(T_o(ix1+1:ix2+1,jy1:jy2,1)+T_o(ix1-1:ix2-1,jy1:jy2,1)-2*T_o(ix1:ix2,jy1:jy2,1))&
  !+((alfa*dt*inRe)/(Pr*dy*dy))*(T_o(ix1:ix2,jy1+1:jy2+1,1)+T_o(ix1:ix2,jy1-1:jy2-1,1)-2*T_o(ix1:ix2,jy1:jy2,1))&
  !-((alfa*dt*u(ix1+1:ix2+1,jy1:jy2,1))/dx)*(T_o(ix1+1:ix2+1,jy1:jy2,1)-T_o(ix1:ix2,jy1:jy2,1))&
  !-((alfa*dt*v(ix1:ix2,jy1+1:jy2+1,1))/dy)*(T_o(ix1:ix2,jy1+1:jy2+1,1)-T_o(ix1:ix2,jy1:jy2,1))



  ! without Loop Vectorization !

  do j=jy1,jy2
     do i=ix1,ix2

     u_conv = (u(i+2,j,1)+u(i+1,j,1))/2.
     v_conv = (v(i,j+2,1)+v(i,j+1,1))/2.

!     u_plus1 = (u(i+2,j,1)+u(i+1,j,1))/2
!     u_mins1 = (u(i,j,1)+u(i-1,j,1))/2

!     v_plus1 = (v(i,j+2,1)+v(i,j+1,1))/2
!     v_mins1 = (v(i,j,1)+v(i,j-1,1))/2

!     if (u_plus1 < u_conv) u_plus1 = 0.
!     if (u_mins1 < u_conv) u_mins1 = 0.
!     if (v_plus1 < v_conv) v_plus1 = 0.
!     if (v_mins1 < v_conv) v_mins1 = 0.

     u_plus = max(u_conv, 0.)
     u_mins = min(u_conv, 0.)

     v_plus = max(v_conv, 0.)
     v_mins = min(v_conv, 0.)

     Tx_plus = T_o(i+1,j,1)-T_o(i,j,1)
     Tx_mins = T_o(i,j,1)-T_o(i-1,j,1)
 !    Tx = T_o(i,j,1)

     Ty_plus = T_o(i,j+1,1)-T_o(i,j,1)
     Ty_mins = T_o(i,j,1)-T_o(i,j-1,1)
 !    Ty = T_o(i,j,1)
     

     T_p(i,j,1) = T_o(i,j,1) + ((alfa*dt*inRe)/(Pr*dx*dx))*(T_o(i+1,j,1)+T_o(i-1,j,1)-2.*T_o(i,j,1))&
                             + ((alfa*dt*inRe)/(Pr*dy*dy))*(T_o(i,j+1,1)+T_o(i,j-1,1)-2.*T_o(i,j,1))&
                             - ((alfa*dt)/dx) * (u_plus*Tx_plus + u_mins*Tx_mins)&
                             - ((alfa*dt)/dy) * (v_plus*Ty_plus + v_mins*Ty_mins)

!    T_p(i,j,1) = T_o(i,j,1) - ((alfa*dt)/dx) * (u_plus*Tx_plus + u_mins*Tx_mins)&
!                            - ((alfa*dt)/dy) * (v_plus*Ty_plus + v_mins*Ty_mins)
     



     end do
  end do 


  do i = ix1,ix2
     T_res = T_res + sum((T_o(i,:,1)-T_p(i,:,1))**2)
 end do

  T_res = sqrt(T_res/size(T_o))

  print *,T_res

end subroutine Heat_Solve
