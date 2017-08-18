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

  real :: u_conv, v_conv
  real :: Tx_plus, Tx_mins, Ty_plus, Ty_mins, Tij
  real :: Txx, Tyy
  real :: coeff, dxm, dxp, dym, dyp, Txp, Txm, Typ, Tym
  real :: tol
  real :: rhoc,rhoxm,rhoym,rhoxp,rhoyp,mdotxm,mdotym,nxm,nym,mdotxp,mdotyp,nxp,nyp

  tol = 0.01

  k = 1

  do j=jy1,jy2
     do i=ix1,ix2

     rhoxm = (smrh(i,j,k) + smrh(i-1,j,k))/2.0d0 - smrh(i,j,k)
     rhoym = (smrh(i,j,k) + smrh(i,j-1,k))/2.0d0 - smrh(i,j,k)

     rhoxp = (smrh(i,j,k) + smrh(i+1,j,k))/2.0d0 - smrh(i,j,k)
     rhoyp = (smrh(i,j,k) + smrh(i,j+1,k))/2.0d0 - smrh(i,j,k)

     mdotxm = mdot(i,j,k)
     mdotym = mdot(i,j,k)
    
     mdotxp = mdot(i,j,k)
     mdotyp = mdot(i,j,k)

     nxm    = nrmx(i,j,k)
     nym    = nrmy(i,j,k)
   
     nxp    = nrmx(i,j,k)
     nyp    = nrmy(i,j,k)

     u_conv = (u(i+1,j,k)+u(i,j,k)+(mdotxm*nxm*rhoxm)+(mdotxp*nxp*rhoxp))/2.
     v_conv = (v(i,j+1,k)+v(i,j,k)+(mdotym*nym*rhoym)+(mdotyp*nyp*rhoyp))/2.

     Tx_plus = T_o(i+1,j,k)
     Tx_mins = T_o(i-1,j,k)

     Ty_plus = T_o(i,j+1,k) 
     Ty_mins = T_o(i,j-1,k)

     Tij = T_o(i,j,k)

     Txp = Tx_plus
     Typ = Ty_plus
     Txm = Tx_mins
     Tym = Ty_mins

     dxp = dx
     dyp = dy
     dxm = dx
     dym = dy

     coeff = inRe/ht_Pr

     ! Case 1 !
     if(s(i,j,k)*s(i+1,j,k) .le. 0.d0) then
     
       th = max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k))))

       Tx_plus = (ht_Tsat-Tij)/th + Tij
       
       Txp = ht_Tsat
       dxp = th*dx

     end if
     ! End of Case 1 !


     ! Case 2 !
     if(s(i,j,k)*s(i-1,j,k).le.0.d0) then
    
       th = max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i-1,j,k))))

       Tx_mins = (ht_Tsat-Tij)/th + Tij

       Txm = ht_Tsat
       dxm = th*dx

     end if
     ! End of Case 2 !

     ! Case 3 !
     if(s(i,j,k)*s(i,j+1,k).le.0.d0) then

       th = max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k))))

       Ty_plus = (ht_Tsat-Tij)/th + Tij

       Typ = ht_Tsat
       dyp = th*dy
 
     end if
     ! End of Case 3 !
 
     ! Case 4 !
     if(s(i,j,k)*s(i,j-1,k).le.0.d0) then

      th = max(tol,abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j-1,k))))

      Ty_mins = (ht_Tsat-Tij)/th + Tij 
       
      Tym = ht_Tsat
      dym = th*dy

     end if
     ! End of Case 4 !  
!_______________________________RHS TERM______________________________________!

    Txx = alph(i,j,k)*(coeff*(Tx_plus-Tij)/dx - coeff*(Tij-Tx_mins)/dx)/dx
    Tyy = alph(i,j,k)*(coeff*(Ty_plus-Tij)/dy - coeff*(Tij-Ty_mins)/dy)/dy

    T_rhs(i,j,k) = - u_conv*(Txp-Txm)/(dxp+dxm) &
                   - v_conv*(Typ-Tym)/(dyp+dym) &
                   + Txx + Tyy

    end do
  end do 

end subroutine Heat_RHS_central
