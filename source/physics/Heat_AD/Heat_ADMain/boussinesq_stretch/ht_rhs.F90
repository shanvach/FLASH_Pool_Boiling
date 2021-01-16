subroutine ht_rhs2d(T, Trhs, u, v, dx, dy, ix1, ix2, jy1, jy2)

  use Heat_AD_data

#include "Heat_AD.h"

  implicit none

#include "constants.h"

  real, dimension(:,:,:), intent(in)  :: T
  real, dimension(:,:,:), intent(out) :: Trhs
  real, dimension(:,:,:), intent(in)  :: u, v
  real, dimension(:,:),   intent(in)  :: dx, dy 
  integer, intent(in) :: ix1, ix2, jy1, jy2

  integer :: i, j

  real :: u_plus, u_mins, v_plus, v_mins
  real :: u_conv, v_conv
  real :: Txplus, Txmins, Typlus, Tymins
  real :: dTdxp,  dTdxm,  dTdyp,  dTdym

#ifdef HEAT_CENTRAL

  do j = jy1, jy2
    do i = ix1, ix2
      
      Txplus = ( T(i+1,j,1) + T(i,j,1)   ) / 2.
      Txmins = ( T(i,j,1)   + T(i-1,j,1) ) / 2.

      Typlus = ( T(i,j+1,1) + T(i,j,1)   ) / 2.
      Tymins = ( T(i,j,1)   + T(i,j-1,1) ) / 2.

      u_conv = u(i+1,j,1) * Txplus - u(i,j,1) * Txmins
      v_conv = v(i,j+1,1) * Typlus - v(i,j,1) * Tymins

      dTdxp = ( T(i+1,j,1) - T(i,j,1)   ) * dx(i,RIGHT_EDGE)
      dTdxm = ( T(i,j,1)   - T(i-1,j,1) ) * dx(i,LEFT_EDGE)
      dTdyp = ( T(i,j+1,1) - T(i,j,1)   ) * dy(j,RIGHT_EDGE)
      dTdym = ( T(i,j,1)   - T(i,j-1,1) ) * dy(j,LEFT_EDGE)
    
      Trhs(i,j,1) =                                       &
                  - u_conv * dx(i,CENTER)                 &
                  - v_conv * dy(j,CENTER)                 &
                  + ht_invsqrtRaPr * (                    &
                       + ht_source_term                   &
                       + (dTdxp - dTdxm) * dx(i,CENTER)   &
                       + (dTdyp - dTdym) * dy(j,CENTER) ) 

    end do
  end do

#endif

#ifdef HEAT_UPWIND

  do j = jy1, jy2
    do i = ix1, ix2
        
      u_conv = ( u(i+1,j,1) + u(i,j,1) ) / 2.
      v_conv = ( v(i,j+1,1) + v(i,j,1) ) / 2.

      u_plus = max(u_conv, 0.)
      u_mins = min(u_conv, 0.)

      v_plus = max(v_conv, 0.)
      v_mins = min(v_conv, 0.)

      Txplus = ( T(i+1,j,1) - T(i,j,1)   ) * dx(i,RIGHT_EDGE)
      Txmins = ( T(i,j,1)   - T(i-1,j,1) ) * dx(i,LEFT_EDGE)

      Typlus = ( T(i,j+1,1) - T(i,j,1)   ) * dy(j,RIGHT_EDGE)
      Tymins = ( T(i,j,1)   - T(i,j-1,1) ) * dy(j,LEFT_EDGE)

      dTdxp = ( T(i+1,j,1) - T(i,j,1)   ) * dx(i,RIGHT_EDGE)
      dTdxm = ( T(i,j,1)   - T(i-1,j,1) ) * dx(i,LEFT_EDGE)
      dTdyp = ( T(i,j+1,1) - T(i,j,1)   ) * dy(j,RIGHT_EDGE)
      dTdym = ( T(i,j,1)   - T(i,j-1,1) ) * dy(j,LEFT_EDGE)
    
      Trhs(i,j,1) =                                       &
                  - ( u_plus * Txmins + u_mins * Txplus ) &                          
                  - ( v_plus * Tymins + v_mins * Typlus ) &                          
                  + ht_invsqrtRaPr * (                    &
                       + ht_source_term                   &
                       + (dTdxp - dTdxm) * dx(i,CENTER)   &
                       + (dTdyp - dTdym) * dy(j,CENTER) ) 

     end do
  end do 

#endif

end subroutine ht_rhs2d


subroutine ht_rhs3d(T, Trhs, u, v, w, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2)

  use Heat_AD_data

#include "Heat_AD.h"

  implicit none
  real, dimension(:,:,:), intent(in)  :: T
  real, dimension(:,:,:), intent(out) :: Trhs
  real, dimension(:,:,:), intent(in)  :: u, v, w
  real, dimension(:,:),   intent(in)  :: dx, dy, dz 
  integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2

  integer :: i, j, k

  real :: u_plus, u_mins, v_plus, v_mins, w_plus, w_mins 
  real :: u_conv, v_conv, w_conv
  real :: Txplus, Txmins, Typlus, Tymins, Tzplus, Tzmins
  real :: dTdxp,  dTdxm,  dTdyp,  dTdym,  dTdzp,  dTdzm

#ifdef HEAT_CENTRAL

  do k = kz1, kz2
    do j = jy1, jy2
      do i = ix1, ix2

        Txplus = ( T(i+1,j,k) + T(i,j,k)   ) / 2.
        Txmins = ( T(i,j,k)   + T(i-1,j,k) ) / 2.

        Typlus = ( T(i,j+1,k) + T(i,j,k)   ) / 2.
        Tymins = ( T(i,j,k)   + T(i,j-1,k) ) / 2.
        
        Tzplus = ( T(i,j,k+1) + T(i,j,k)   ) / 2.
        Tzmins = ( T(i,j,k)   + T(i,j,k-1) ) / 2.

        u_conv = u(i+1,j,k) * Txplus - u(i,j,k) * Txmins
        v_conv = v(i,j+1,k) * Typlus - v(i,j,k) * Tymins
        w_conv = w(i,j,k+1) * Tzplus - w(i,j,k) * Tzmins

        dTdxp = ( T(i+1,j,k) - T(i,j,k)   ) * dx(i,RIGHT_EDGE)
        dTdxm = ( T(i,j,k)   - T(i-1,j,k) ) * dx(i,LEFT_EDGE)
        dTdyp = ( T(i,j+1,k) - T(i,j,k)   ) * dy(j,RIGHT_EDGE)
        dTdym = ( T(i,j,k)   - T(i,j-1,k) ) * dy(j,LEFT_EDGE)
        dTdzp = ( T(i,j,k+1) - T(i,j,k)   ) * dz(k,RIGHT_EDGE)
        dTdzm = ( T(i,j,k)   - T(i,j,k-1) ) * dz(k,LEFT_EDGE)

        Trhs(i,j,k) =                                       &
                    - u_conv * dx(i,CENTER)                 &
                    - v_conv * dy(j,CENTER)                 &
                    - w_conv * dz(k,CENTER)                 &
                    + ht_invsqrtRaPr * (                    &
                         + ht_source_term                   &
                         + (dTdxp - dTdxm) * dx(i,CENTER)   &
                         + (dTdyp - dTdym) * dy(j,CENTER)   &
                         + (dTdzp - dTdzm) * dz(k,CENTER) )

      end do
    end do
  end do

#endif


#ifdef HEAT_UPWIND

  do k=kz1, kz2
    do j=jy1, jy2
      do i=ix1, ix2

        u_conv = ( u(i+1,j,k) + u(i,j,k) ) / 2.
        v_conv = ( v(i,j+1,k) + v(i,j,k) ) / 2.
        w_conv = ( w(i,j,k+1) + w(i,j,k) ) / 2.

        u_plus = max(u_conv, 0.)
        u_mins = min(u_conv, 0.)

        v_plus = max(v_conv, 0.)
        v_mins = min(v_conv, 0.)

        w_plus = max(w_conv, 0.)
        w_mins = min(w_conv, 0.)

        Txplus = ( T(i+1,j,k) - T(i,j,k)   ) * dx(i,RIGHT_EDGE)
        Txmins = ( T(i,j,k)   - T(i-1,j,k) ) * dx(i,LEFT_EDGE)

        Typlus = ( T(i,j+1,k) - T(i,j,k)   ) * dy(j,RIGHT_EDGE)
        Tymins = ( T(i,j,k)   - T(i,j-1,k) ) * dy(j,LEFT_EDGE)

        Tzplus = ( T(i,j,k+1) - T(i,j,k)   ) * dz(k,RIGHT_EDGE)
        Tzmins = ( T(i,j,k)   - T(i,j,k-1) ) * dz(k,LEFT_EDGE)

        dTdxp = ( T(i+1,j,k) - T(i,j,k)   ) * dx(i,RIGHT_EDGE)
        dTdxm = ( T(i,j,k)   - T(i-1,j,k) ) * dx(i,LEFT_EDGE)
        dTdyp = ( T(i,j+1,k) - T(i,j,k)   ) * dy(j,RIGHT_EDGE)
        dTdym = ( T(i,j,k)   - T(i,j-1,k) ) * dy(j,LEFT_EDGE)
        dTdzp = ( T(i,j,k+1) - T(i,j,k)   ) * dz(k,RIGHT_EDGE)
        dTdzm = ( T(i,j,k)   - T(i,j,k-1) ) * dz(k,LEFT_EDGE)
    
        Trhs(i,j,k) =                                       &
                    - ( u_plus * Txmins + u_mins * Txplus)  &                          
                    - ( v_plus * Tymins + v_mins * Typlus)  &                          
                    - ( w_plus * Tzmins + w_mins * Tzplus)  &                          
                    + ht_invsqrtRaPr * (                    &
                         + ht_source_term                   &
                         + (dTdxp - dTdxm) * dx(i,CENTER)   &
                         + (dTdyp - dTdym) * dy(j,CENTER)   & 
                         + (dTdzp - dTdzm) * dz(k,CENTER) ) 

        end do
     end do
  end do 

#endif

end subroutine ht_rhs3d
