subroutine Heat_RHS_3D(T_rhs, T_o, u, v, w,dx, dy, dz,inRe, ix1,ix2, jy1,jy2,&
                       kz1,kz2,rho1x,rho2x,rho1y,rho2y,rho1z,rho2z,alph,pf,s,mdot,nrmx,nrmy,nrmz,smrh,curv)

  use Heat_AD_data
  use Multiphase_data, only: mph_cp2,mph_thco2, mph_rho2,mph_rho1

#include "Heat_AD.h"

!#define FIRST_ORDER_UPWIND
#define SECOND_ORDER_CENTRAL

  implicit none
  real, dimension(:,:,:), intent(inout) :: T_rhs
  real, dimension(:,:,:), intent(in) :: T_o
  real, dimension(:,:,:), intent(in) :: u,v,w
  real, intent(in) :: dx, dy, dz, inRe
  integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
  real, dimension(:,:,:),intent(in) :: rho1x,rho2x,rho1y,rho2y,alph,rho1z,rho2z
  real, dimension(:,:,:),intent(in) :: pf,s,mdot,nrmx,nrmy,nrmz,smrh,curv

  real :: T_res,th

  integer :: i,j,k

  real :: u_plus, u_mins, v_plus, v_mins, w_plus, w_mins, u_conv, v_conv, w_conv
  real :: Tx_plus, Tx_mins, Ty_plus, Ty_mins, Tz_plus, Tz_mins, Tij
  real :: Txx, Tyy, Tzz
  real :: coeff 
  real :: tol
  real :: rhoc,rhoxm,rhoym,rhozm,rhoxp,rhoyp,rhozp
  real :: ul,ur,vl,vr,wl,wr
  real :: Txl,Txr,Tyl,Tyr,Tzl,Tzr

  tol = 0.01

 do k=kz1,kz2
  do j=jy1,jy2
     do i=ix1,ix2

     rhoxm = (smrh(i,j,k) + smrh(i-1,j,k))/2.0d0 - smrh(i,j,k)
     rhoym = (smrh(i,j,k) + smrh(i,j-1,k))/2.0d0 - smrh(i,j,k)
     rhozm = (smrh(i,j,k) + smrh(i,j,k-1))/2.0d0 - smrh(i,j,k)

     rhoxp = (smrh(i,j,k) + smrh(i+1,j,k))/2.0d0 - smrh(i,j,k)
     rhoyp = (smrh(i,j,k) + smrh(i,j+1,k))/2.0d0 - smrh(i,j,k)
     rhozp = (smrh(i,j,k) + smrh(i,j,k+1))/2.0d0 - smrh(i,j,k)

     u_conv = (u(i+1,j,k)+u(i,j,k)+(mdot(i,j,k)*nrmx(i,j,k)*rhoxm)+(mdot(i,j,k)*nrmx(i,j,k)*rhoxp))/2.
     v_conv = (v(i,j+1,k)+v(i,j,k)+(mdot(i,j,k)*nrmy(i,j,k)*rhoym)+(mdot(i,j,k)*nrmy(i,j,k)*rhoyp))/2.
     w_conv = (w(i,j,k+1)+w(i,j,k)+(mdot(i,j,k)*nrmz(i,j,k)*rhozm)+(mdot(i,j,k)*nrmz(i,j,k)*rhozp))/2.

     u_plus = max(u_conv, 0.)
     u_mins = min(u_conv, 0.)

     v_plus = max(v_conv, 0.)
     v_mins = min(v_conv, 0.)   

     w_plus = max(w_conv, 0.)
     w_mins = min(w_conv, 0.)

     ul = u(i,j,k)   + mdot(i,j,k)*nrmx(i,j,k)*rhoxm
     ur = u(i+1,j,k) + mdot(i,j,k)*nrmx(i,j,k)*rhoxp
     vl = v(i,j,k)   + mdot(i,j,k)*nrmy(i,j,k)*rhoym
     vr = v(i,j+1,k) + mdot(i,j,k)*nrmy(i,j,k)*rhoyp
     wl = w(i,j,k)   + mdot(i,j,k)*nrmz(i,j,k)*rhozm
     wr = w(i,j,k+1) + mdot(i,j,k)*nrmz(i,j,k)*rhozp

     Tx_plus = T_o(i+1,j,k)
     Tx_mins = T_o(i-1,j,k)

     Ty_plus = T_o(i,j+1,k) 
     Ty_mins = T_o(i,j-1,k)

     Tz_plus = T_o(i,j,k+1)
     Tz_mins = T_o(i,j,k-1)

     Tij = T_o(i,j,k)

     coeff = inRe/ht_Pr

     ! Case 1 !
     if(s(i,j,k)*s(i+1,j,k) .le. 0.d0) then

       th = max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k))))

       if(s(i,j,k)*s(i-1,j,k) .le. 0.d0) then
       Tx_plus = (ht_Tsat-Tij)/th + Tij

       else
       Tx_plus = (2*ht_Tsat + (2*th*th - 2)*Tij + (-th*th + th)*T_o(i-1,j,k))/(th + th*th)

       end if

     end if
     ! End of Case 1 !


     ! Case 2 !
     if(s(i,j,k)*s(i-1,j,k) .le. 0.d0) then

       th = max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i-1,j,k))))

       if(s(i,j,k)*s(i+1,j,k) .le. 0.d0) then
       Tx_mins = (ht_Tsat-Tij)/th + Tij

       else
       Tx_mins = (2*ht_Tsat + (2*th*th - 2)*Tij + (-th*th + th)*T_o(i+1,j,k))/(th + th*th)

       end if

     end if
     ! End of Case 2 !

     ! Case 3 !
     if(s(i,j,k)*s(i,j+1,k) .le. 0.d0) then

       th = max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k))))

       if(s(i,j,k)*s(i,j-1,k) .le. 0.0d0) then
       Ty_plus = (ht_Tsat-Tij)/th + Tij

       else
       Ty_plus = (2*ht_Tsat + (2*th*th - 2)*Tij + (-th*th + th)*T_o(i,j-1,k))/(th + th*th)

       end if

     end if
     ! End of Case 3 !

     ! Case 4 !
     if(s(i,j,k)*s(i,j-1,k) .le. 0.d0) then

      th = max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j-1,k))))

      if(s(i,j,k)*s(i,j+1,k) .le. 0.d0) then
      Ty_mins = (ht_Tsat-Tij)/th + Tij

      else
      Ty_mins = (2*ht_Tsat + (2*th*th - 2)*Tij + (-th*th + th)*T_o(i,j+1,k))/(th + th*th)

      end if

     end if
     ! End of Case 4 ! 

     ! Case 5 !
     if(s(i,j,k)*s(i,j,k+1) .le. 0.d0) then

      th = max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j,k+1))))

      if(s(i,j,k)*s(i,j,k-1) .le. 0.d0) then
      Tz_plus = (ht_Tsat-Tij)/th + Tij

      else
      Tz_plus = (2*ht_Tsat + (2*th*th - 2)*Tij + (-th*th + th)*T_o(i,j,k-1))/(th + th*th)

      end if

     end if
     ! End of Case 5 !
 
     ! Case 6 !
     if(s(i,j,k)*s(i,j,k-1) .le. 0.d0) then

      th = max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j,k-1))))

      if(s(i,j,k)*s(i,j,k+1) .le. 0.d0) then
      Tz_mins = (ht_Tsat-Tij)/th + Tij

      else
      Tz_mins = (2*ht_Tsat + (2*th*th - 2)*Tij + (-th*th + th)*T_o(i,j,k+1))/(th + th*th)

      end if

     end if
     ! End of Case 6 !

!_________________________________RHS TERM____________________________________!

    Txx = alph(i,j,k)*(coeff*(Tx_plus-Tij)/dx - coeff*(Tij-Tx_mins)/dx)/dx
    Tyy = alph(i,j,k)*(coeff*(Ty_plus-Tij)/dy - coeff*(Tij-Ty_mins)/dy)/dy
    Tzz = alph(i,j,k)*(coeff*(Tz_plus-Tij)/dz - coeff*(Tij-Tz_mins)/dz)/dz

    Txl  = (Tx_mins + Tij)/2.0d0
    Txr  = (Tx_plus + Tij)/2.0d0
    Tyl  = (Ty_mins + Tij)/2.0d0
    Tyr  = (Ty_plus + Tij)/2.0d0
    Tzl  = (Tz_mins + Tij)/2.0d0
    Tzr  = (Tz_plus + Tij)/2.0d0

#ifdef FIRST_ORDER_UPWIND
    T_rhs(i,j,k) = ((-(u_plus*(Tij-Tx_mins)/dx+u_mins*(Tx_plus-Tij)/dx)&
                     -(v_plus*(Tij-Ty_mins)/dy+v_mins*(Ty_plus-Tij)/dy)&
                     -(w_plus*(Tij-Tz_mins)/dz+w_mins*(Tz_plus-Tij)/dz))&
                     +(Txx+Tyy+Tzz))      
#endif

#ifdef SECOND_ORDER_CENTRAL
    !T_rhs(i,j,k) = - u_conv*(Tx_plus-Tx_mins)/(2*dx) &
    !               - v_conv*(Ty_plus-Ty_mins)/(2*dy) &
    !               - w_conv*(Tz_plus-Tz_mins)/(2*dz) &
    !               + Txx + Tyy + Tzz

    T_rhs(i,j,k) =  - (ur*Txr-ul*Txl)/dx &
                    - (vr*Tyr-vl*Tyl)/dy &
                    - (wr*Tzr-wl*Tzl)/dz &
                    + Txx + Tyy + Tzz
#endif

     end do
  end do 
end do

end subroutine Heat_RHS_3D
