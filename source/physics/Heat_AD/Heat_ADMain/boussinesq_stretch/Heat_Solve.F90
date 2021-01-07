subroutine Heat_Solve2d(T_p, T_o, u, v, dt, dx, dy, ix1,ix2, jy1, jy2)

  use Heat_AD_data

#include "Heat_AD.h"

  implicit none

#include "constants.h"

  real, dimension(:,:,:), intent(inout) :: T_p
  real, dimension(:,:,:), intent(in) :: T_o
  real, dimension(:,:,:), intent(in) :: u,v
  real, intent(in) :: dt
  real, dimension(:,:),   intent(in) :: dx, dy 
  integer, intent(in) :: ix1, ix2, jy1, jy2

  integer :: i,j

  real :: u_plus, u_mins, v_plus, v_mins
  real :: u_conv, v_conv
  real :: Txplus, Txmins, Typlus, Tymins
  real :: dTdxp,  dTdxm,  dTdyp,  dTdym

#ifdef HEAT_CENTRAL

  ! without Loop Vectorization ! (Second Order Central)
  do j = jy1,jy2
    do i = ix1,ix2
      
      u_conv = (u(i+1,j,1) + u(i,j,1))/2.
      v_conv = (v(i,j+1,1) + v(i,j,1))/2.

      Txplus = (T_o(i+1,j,1) + T_o(i,j,1))/2.
      Txmins = (T_o(i,j,1) + T_o(i-1,j,1))/2.

      Typlus = (T_o(i,j+1,1) + T_o(i,j,1))/2.
      Tymins = (T_o(i,j,1) + T_o(i,j-1,1))/2.

      dTdxp = (T_o(i+1,j,1) - T_o(i,j,1)) * dx(i,RIGHT_EDGE)
      dTdxm = (T_o(i,j,1) - T_o(i-1,j,1)) * dx(i,LEFT_EDGE)
      dTdyp = (T_o(i,j+1,1) - T_o(i,j,1)) * dy(j,RIGHT_EDGE)
      dTdym = (T_o(i,j,1) - T_o(i,j-1,1)) * dy(j,LEFT_EDGE)
    
      T_p(i,j,1) =                                                  &
                  + T_o(i,j,1)                                      &
                  - dt * u_conv * (Txplus - Txmins) * dx(i,CENTER)  &                          
                  - dt * v_conv * (Typlus - Tymins) * dy(j,CENTER)  &                          
                  + dt * ht_invsqrtRaPr * (                         &
                       + ht_source_term                             &
                       + (dTdxp - dTdxm) * dx(i,CENTER)             &
                       + (dTdyp - dTdym) * dy(j,CENTER)) 

    end do
  end do

#endif


#ifdef HEAT_UPWIND
  ! without Loop Vectorization ! (First Order Upwind)

  do j=jy1,jy2
    do i=ix1,ix2
        
      u_conv = (u(i+1,j,1) + u(i,j,1))/2.
      v_conv = (v(i,j+1,1) + v(i,j,1))/2.

      u_plus = max(u_conv, 0.)
      u_mins = min(u_conv, 0.)

      v_plus = max(v_conv, 0.)
      v_mins = min(v_conv, 0.)

      Txplus = (T_o(i+1,j,1) - T_o(i,j,1)) * dx(i,RIGHT_EDGE)
      Txmins = (T_o(i,j,1) - T_o(i-1,j,1)) * dx(i,LEFT_EDGE)

      Typlus = (T_o(i,j+1,1) - T_o(i,j,1)) * dy(j,RIGHT_EDGE)
      Tymins = (T_o(i,j,1) - T_o(i,j-1,1)) * dy(j,LEFT_EDGE)

      dTdxp = (T_o(i+1,j,1) - T_o(i,j,1)) * dx(i,RIGHT_EDGE)
      dTdxm = (T_o(i,j,1) - T_o(i-1,j,1)) * dx(i,LEFT_EDGE)
      dTdyp = (T_o(i,j+1,1) - T_o(i,j,1)) * dy(j,RIGHT_EDGE)
      dTdym = (T_o(i,j,1) - T_o(i,j-1,1)) * dy(j,LEFT_EDGE)
    
      T_p(i,j,1) =                                                  &
                  + T_o(i,j,1)                                      &
                  - dt * (u_plus*Txmins + u_mins*Txplus)            &                          
                  - dt * (v_plus*Tymins + v_mins*Typlus)            &                          
                  + dt * ht_invsqrtRaPr * (                         &
                       + ht_source_term                             &
                       + (dTdxp - dTdxm) * dx(i,CENTER)             &
                       + (dTdyp - dTdym) * dy(j,CENTER)) 

     end do
  end do 

#endif

end subroutine Heat_Solve2d


subroutine Heat_Solve3d(T_p, T_o, u, v, w, dt, dx, dy, dz, &
                        ix1,ix2, jy1, jy2, kz1, kz2)

  use Heat_AD_data

#include "Heat_AD.h"

  implicit none
  real, dimension(:,:,:), intent(inout) :: T_p
  real, dimension(:,:,:), intent(in) :: T_o
  real, dimension(:,:,:), intent(in) :: u,v,w
  real, intent(in) :: dt
  real, dimension(:,:),   intent(in) :: dx, dy, dz 
  integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2

  integer :: i,j,k

  real :: u_plus, u_mins, v_plus, v_mins, w_plus, w_mins 
  real :: u_conv, v_conv, w_conv
  real :: Txplus, Txmins, Typlus, Tymins, Tzplus, Tzmins
  real :: dTdxp,  dTdxm,  dTdyp,  dTdym,  dTdzp,  dTdzm

#ifdef HEAT_CENTRAL

  ! without Loop Vectorization ! (Second Order Central)
  do k = kz1,kz2
    do j = jy1,jy2
      do i = ix1,ix2

        u_conv = (u(i+1,j,k) + u(i,j,k))/2.
        v_conv = (v(i,j+1,k) + v(i,j,k))/2.
        w_conv = (w(i,j,k+1) + w(i,j,k))/2.

        Txplus = (T_o(i+1,j,k) + T_o(i,j,k))/2.
        Txmins = (T_o(i,j,k) + T_o(i-1,j,k))/2.

        Typlus = (T_o(i,j+1,k) + T_o(i,j,k))/2.
        Tymins = (T_o(i,j,k) + T_o(i,j-1,k))/2.

        Tzplus = (T_o(i,j,k+1) + T_o(i,j,k))/2.
        Tzmins = (T_o(i,j,k) + T_o(i,j,k-1))/2.

        dTdxp = (T_o(i+1,j,k) - T_o(i,j,k)) * dx(i,RIGHT_EDGE)
        dTdxm = (T_o(i,j,k) - T_o(i-1,j,k)) * dx(i,LEFT_EDGE)
        dTdyp = (T_o(i,j+1,k) - T_o(i,j,k)) * dy(j,RIGHT_EDGE)
        dTdym = (T_o(i,j,k) - T_o(i,j-1,k)) * dy(j,LEFT_EDGE)
        dTdzp = (T_o(i,j,k+1) - T_o(i,j,k)) * dz(k,RIGHT_EDGE)
        dTdzm = (T_o(i,j,k) - T_o(i,j,k-1)) * dz(k,LEFT_EDGE)

        T_p(i,j,k) =                                                  &
                    + T_o(i,j,k)                                      &
                    - dt * u_conv * (Txplus - Txmins) * dx(i,CENTER)  &
                    - dt * v_conv * (Typlus - Tymins) * dy(j,CENTER)  &
                    - dt * w_conv * (Tzplus - Tzmins) * dz(k,CENTER)  &
                    + dt * ht_invsqrtRaPr * (                         &
                         + ht_source_term                             &
                         + (dTdxp - dTdxm) * dx(i,CENTER)             &
                         + (dTdyp - dTdym) * dy(j,CENTER)             &
                         + (dTdzp - dTdzm) * dz(k,CENTER))

      end do
    end do
  end do

#endif


#ifdef HEAT_UPWIND

  ! without Loop Vectorization ! (First Order Upwind)
  do k=kz1,kz2
    do j=jy1,jy2
      do i=ix1,ix2

        u_conv = (u(i+1,j,k) + u(i,j,k))/2.
        v_conv = (v(i,j+1,k) + v(i,j,k))/2.
        w_conv = (w(i,j,k+1) + w(i,j,k))/2.

        u_plus = max(u_conv, 0.)
        u_mins = min(u_conv, 0.)

        v_plus = max(v_conv, 0.)
        v_mins = min(v_conv, 0.)

        w_plus = max(w_conv, 0.)
        w_mins = min(w_conv, 0.)

        Tx_plus = (T_o(i+1,j,k) - T_o(i,j,k)) * dx(i,RIGHT_EDGE)
        Tx_mins = (T_o(i,j,k) - T_o(i-1,j,k)) * dx(i,LEFT_EGDE)

        Ty_plus = (T_o(i,j+1,k) - T_o(i,j,k)) * dy(j,RIGHT_EDGE)
        Ty_mins = (T_o(i,j,k) - T_o(i,j-1,k)) * dy(j,LEFT_EDGE)

        Tz_plus = (T_o(i,j,k+1) - T_o(i,j,k)) * dz(k,RIGHT_EDGE)
        Tz_mins = (T_o(i,j,k) - T_o(i,j,k-1)) * dz(k,LEFT_EDGE)

        dTdxp = (T_o(i+1,j,k) - T_o(i,j,k)) * dx(i,RIGHT_EDGE)
        dTdxm = (T_o(i,j,k) - T_o(i-1,j,k)) * dx(i,LEFT_EDGE)
        dTdyp = (T_o(i,j+1,k) - T_o(i,j,k)) * dy(j,RIGHT_EDGE)
        dTdym = (T_o(i,j,k) - T_o(i,j-1,k)) * dy(j,LEFT_EDGE)
        dTdzp = (T_o(i,j,k+1) - T_o(i,j,k)) * dz(k,RIGHT_EDGE)
        dTdzm = (T_o(i,j,k) - T_o(i,j,k-1)) * dz(k,LEFT_EDGE)
    
        T_p(i,j,k) =                                                  &
                    + T_o(i,j,k)                                      &
                    - dt * (u_plus*Txmins + u_mins*Txplus)            &                          
                    - dt * (v_plus*Tymins + v_mins*Typlus)            &                          
                    - dt * (w_plus*Twmins + w_mins*Twplus)            &                          
                    + dt * ht_invsqrtRaPr * (                         &
                         + ht_source_term                             &
                         + (dTdxp - dTdxm) * dx(i,CENTER)             &
                         + (dTdyp - dTdym) * dy(j,CENTER)             & 
                         + (dTdwp - dTdwm) * dz(k,CENTER)) 

        end do
     end do
  end do 

#endif

end subroutine Heat_Solve3d
