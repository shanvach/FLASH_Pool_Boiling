subroutine Heat_RHS_central(T_rhs, T_o, u, v, dx, dy, dz,inRe, ix1,ix2, jy1,jy2,&
                    rho1x,rho2x,rho1y,rho2y,alph,pf,s,mdot,nrmx,nrmy,smrh,curv)

  use Heat_AD_data
  use Multiphase_data, only: mph_cp2,mph_thco2, mph_rho2,mph_rho1

#include "Heat_AD.h"

  implicit none
  real, dimension(:,:,:), intent(inout) :: T_rhs
  real, dimension(:,:,:), intent(in) :: T_o
  real, dimension(:,:,:), intent(in) :: u,v
  real, intent(in) :: dx, dy, dz, inRe
  integer, intent(in) :: ix1, ix2, jy1, jy2
  real, dimension(:,:,:),intent(in) :: rho1x,rho2x,rho1y,rho2y,alph
  real, dimension(:,:,:),intent(in) :: pf,s,mdot,nrmx,nrmy,smrh,curv

  real :: T_res,th

  integer :: i,j,k

  real :: u_conv, v_conv, ul, ur, vl, vr
  real :: Tx_plus, Tx_mins, Ty_plus, Ty_mins, Tij
  real :: Txx, Tyy, Txl, Txr, Tyl, Tyr
  real :: coeff
  real :: tol
  real :: rhoc,rhoxm,rhoym,rhoxp,rhoyp

  tol = 0.01

  k = 1

  do j=jy1,jy2
     do i=ix1,ix2

     rhoxm = (smrh(i,j,k) + smrh(i-1,j,k))/2.0d0 - smrh(i,j,k)
     rhoym = (smrh(i,j,k) + smrh(i,j-1,k))/2.0d0 - smrh(i,j,k)

     rhoxp = (smrh(i,j,k) + smrh(i+1,j,k))/2.0d0 - smrh(i,j,k)
     rhoyp = (smrh(i,j,k) + smrh(i,j+1,k))/2.0d0 - smrh(i,j,k)

     u_conv = (u(i+1,j,k)+u(i,j,k)+(mdot(i,j,k)*nrmx(i,j,k)*rhoxm)+(mdot(i,j,k)*nrmx(i,j,k)*rhoxp))/2.
     v_conv = (v(i,j+1,k)+v(i,j,k)+(mdot(i,j,k)*nrmy(i,j,k)*rhoym)+(mdot(i,j,k)*nrmy(i,j,k)*rhoyp))/2.

     ul = u(i,j,k)   + mdot(i,j,k)*nrmx(i,j,k)*rhoxm
     ur = u(i+1,j,k) + mdot(i,j,k)*nrmx(i,j,k)*rhoxp
     vl = v(i,j,k)   + mdot(i,j,k)*nrmy(i,j,k)*rhoym
     vr = v(i,j+1,k) + mdot(i,j,k)*nrmy(i,j,k)*rhoyp

     Tx_plus = T_o(i+1,j,k)
     Tx_mins = T_o(i-1,j,k)

     Ty_plus = T_o(i,j+1,k) 
     Ty_mins = T_o(i,j-1,k)

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
!_______________________________RHS TERM______________________________________!

    Txx = alph(i,j,k)*(coeff*(Tx_plus-Tij)/dx - coeff*(Tij-Tx_mins)/dx)/dx
    Tyy = alph(i,j,k)*(coeff*(Ty_plus-Tij)/dy - coeff*(Tij-Ty_mins)/dy)/dy

    Txl  = (Tx_mins + Tij)/2.0d0
    Txr  = (Tx_plus + Tij)/2.0d0
    Tyl  = (Ty_mins + Tij)/2.0d0
    Tyr  = (Ty_plus + Tij)/2.0d0

    !T_rhs(i,j,k) =  - u_conv*(Tx_plus-Tx_mins)/(2*dx) &
    !                - v_conv*(Ty_plus-Ty_mins)/(2*dy) &
    !                + Txx + Tyy

    T_rhs(i,j,k) =  - (ur*Txr-ul*Txl)/dx &
                    - (vr*Tyr-vl*Tyl)/dy &
                    + Txx + Tyy

    end do
  end do 

end subroutine Heat_RHS_central
