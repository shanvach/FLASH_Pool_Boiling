subroutine Heat_RHS(T_rhs, T_o, u, v, dx, dy, dz,inRe, ix1,ix2, jy1,jy2,&
                    rho1x,rho2x,rho1y,rho2y,thco,cp,pf,s,mdot,nrmx,nrmy,smrh)

  use Heat_AD_data
  use Multiphase_data, only: mph_cp2,mph_thco2, mph_rho2,mph_rho1

#include "Heat_AD.h"

  implicit none
  real, dimension(:,:,:), intent(inout) :: T_rhs
  real, dimension(:,:,:), intent(in) :: T_o
  real, dimension(:,:,:), intent(in) :: u,v
  real, intent(in) :: dx, dy, dz, inRe
  integer, intent(in) :: ix1, ix2, jy1, jy2
  real, dimension(:,:,:),intent(in) :: rho1x,rho2x,rho1y,rho2y,thco,cp
  real, dimension(:,:,:),intent(in) :: pf,s,mdot,nrmx,nrmy,smrh

  real :: T_res,Mdensx,Mdensy,th,dxp,dxm,dyp,dym

  integer :: i,j,k

  real :: u_plus, u_mins, v_plus, v_mins, u_conv, v_conv
  real :: Tx_plus, Tx_mins, Ty_plus, Ty_mins, Tij, Tipj, Timj, Tijp, Tijm
  real :: Txx, Tyy
  real :: alphax_plus, alphax_mins, alphay_plus, alphay_mins, alpha_interface
  real :: alpha
  real :: tol
  real :: rhoc,rhoxm,rhoym,rhoxp,rhoyp,mdotxm,mdotym,nxm,nym,mdotxp,mdotyp,nxp,nyp

  tol = 0.01

  k = 1

  do j=jy1,jy2
     do i=ix1,ix2

     !rhoc = (rho1x(i+1,j,k)+rho1x(i,j,k)+&
     !        rho1y(i,j+1,k)+rho2y(i,j,k))/2.0

     rhoxm = rho1x(i,j,k) + rho2x(i,j,k) - smrh(i,j,k)
     rhoym = rho1y(i,j,k) + rho2y(i,j,k) - smrh(i,j,k)

     rhoxp = rho1x(i+1,j,k) + rho2x(i+1,j,k) - smrh(i,j,k)
     rhoyp = rho1y(i,j+1,k) + rho2y(i,j+1,k) - smrh(i,j,k)

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

     u_plus = max(u_conv, 0.)
     u_mins = min(u_conv, 0.)

     v_plus = max(v_conv, 0.)
     v_mins = min(v_conv, 0.)   

     Tx_plus = T_o(i+1,j,k)
     Tx_mins = T_o(i-1,j,k)

     Ty_plus = T_o(i,j+1,k) 
     Ty_mins = T_o(i,j-1,k)

     Tij = T_o(i,j,k)

     ! Case 1 !
     if(s(i,j,k)*s(i+1,j,k) .le. 0.d0) then
     
       if (abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k))) .gt. tol) then

       th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k)))

       else

       th = tol

       end if 

       Tx_plus = (ht_Tsat-Tij)/th + Tij

     end if
     ! End of Case 1 !


     ! Case 2 !
     if(s(i,j,k)*s(i-1,j,k).le.0.d0) then
    
       if(abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i-1,j,k))) .gt. tol) then

       th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i-1,j,k)))

       else

       th = tol
      
       end if

       Tx_mins = (ht_Tsat-Tij)/th + Tij

     end if
     ! End of Case 2 !

    ! Case 3 !
    if(s(i,j,k)*s(i,j+1,k).le.0.d0) then

      if (abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k))) .gt. tol) then

      th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k)))

      else

      th = tol

      end if
  
      Ty_plus = (ht_Tsat-Tij)/th + Tij

   end if
    ! End of Case 3 !
 
    ! Case 4 !
    if(s(i,j,k)*s(i,j-1,k).le.0.d0) then

      if (abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j-1,k))) .gt. tol) then

      th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j-1,k)))

      else

      th = tol

      end if

      Ty_mins = (ht_Tsat-Tij)/th + Tij 

    end if
    ! End of Case 4 !

    alphax_plus = (thco(i,j,k)/cp(i,j,k))*(inRe/ht_Pr)
    alphax_mins = (thco(i,j,k)/cp(i,j,k))*(inRe/ht_Pr)
    alphay_plus = (thco(i,j,k)/cp(i,j,k))*(inRe/ht_Pr)
    alphay_mins = (thco(i,j,k)/cp(i,j,k))*(inRe/ht_Pr)

!_______________________________RHS TERM______________________________________!

    Txx = (alphax_plus*(Tx_plus-Tij)/dx - alphax_mins*(Tij-Tx_mins)/dx)/dx
    Tyy = (alphay_plus*(Ty_plus-Tij)/dy - alphay_mins*(Tij-Ty_mins)/dy)/dy

    T_rhs(i,j,k) = ((-(u_plus*(Tij-Tx_mins)/dx+u_mins*(Tx_plus-Tij)/dx)&
                     -(v_plus*(Tij-Ty_mins)/dy+v_mins*(Ty_plus-Tij)/dy))&
                     +(Txx +Tyy))      

     end do
  end do 

end subroutine Heat_RHS
