subroutine Heat_RHS_3D(T_rhs, T_o, u, v, dx, dy, dz,inRe, ix1,ix2, jy1,jy2,&
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
  real :: rhoxm,rhoym,rhoxp,rhoyp,mdotxm,mdotym,nxm,nym,mdotxp,mdotyp,nxp,nyp

  tol = 0.01

!#ifdef HEAT_CENTRAL
!
!  ! with Loop Vectorization ! (Second Order Central)
!
!  T_rhs(ix1:ix2,jy1:jy2,1) = ((inRe)/(ht_Pr*dx*dx))*(T_o(ix1+1:ix2+1,jy1:jy2,1)+T_o(ix1-1:ix2-1,jy1:jy2,1)-2*T_o(ix1:ix2,jy1:jy2,1))&
!  +((inRe)/(ht_Pr*dy*dy))*(T_o(ix1:ix2,jy1+1:jy2+1,1)+T_o(ix1:ix2,jy1-1:jy2-1,1)-2*T_o(ix1:ix2,jy1:jy2,1))&
!  -(((u(ix1+1:ix2+1,jy1:jy2,1) + u(ix1:ix2,jy1:jy2,1))/2)/(2*dx))*(T_o(ix1+1:ix2+1,jy1:jy2,1)-T_o(ix1-1:ix2-1,jy1:jy2,1))&
!  -(((v(ix1:ix2,jy1+1:jy2+1,1) + v(ix1:ix2,jy1:jy2,1))/2)/(2*dy))*(T_o(ix1:ix2,jy1+1:jy2+1,1)-T_o(ix1:ix2,jy1-1:jy2-1,1))
!
!#endif

!#ifdef HEAT_UPWIND
  ! without Loop Vectorization ! (First Order Upwind)

  k = 1

  do j=jy1,jy2
     do i=ix1,ix2

     rhoxm = 2.d0/(smrh(i-1,j,k)+smrh(i,j,k)) - 1./smrh(i,j,k)
     rhoym = 2.d0/(smrh(i,j-1,k)+smrh(i,j,k)) - 1./smrh(i,j,k)

     rhoxp = 2.d0/(smrh(i+1,j,k)+smrh(i,j,k)) - 1./smrh(i,j,k)
     rhoyp = 2.d0/(smrh(i,j+1,k)+smrh(i,j,k)) - 1./smrh(i,j,k)

     !rhoxm = ((rho1x(i,j,k)+rho2x(i,j,k))/mph_rho2) - 1./smrh(i,j,k)
     !rhoym = ((rho1y(i,j,k)+rho2y(i,j,k))/mph_rho2) - 1./smrh(i,j,k)

     !rhoxp = ((rho1x(i+1,j,k)+rho2x(i+1,j,k))/mph_rho2) - 1./smrh(i,j,k)
     !rhoyp = ((rho1y(i,j+1,k)+rho2y(i,j+1,k))/mph_rho2) - 1./smrh(i,j,k)

     !rhoxm = rhoxm/mph_rho2
     !rhoym = rhoym/mph_rho2
     !rhoxp = rhoxp/mph_rho2
     !rhoyp = rhoyp/mph_rho2

     mdotxm = (mdot(i,j,k)+mdot(i-1,j,k))/2.
     mdotym = (mdot(i,j,k)+mdot(i,j-1,k))/2.

     nxm = (nrmx(i,j,k)+nrmx(i-1,j,k))/2.
     nym = (nrmy(i,j,k)+nrmy(i,j-1,k))/2. 

     mdotxp = (mdot(i,j,k)+mdot(i+1,j,k))/2.
     mdotyp = (mdot(i,j,k)+mdot(i,j+1,k))/2.

     nxp = (nrmx(i,j,k)+nrmx(i+1,j,k))/2.
     nyp = (nrmy(i,j,k)+nrmy(i,j+1,k))/2.

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

     !Tipj = T_o(i+1,j,k)
     !Timj = T_o(i-1,j,k)

     !Tijp = T_o(i,j+1,k)
     !Tijm = T_o(i,j-1,k)

     Tij = T_o(i,j,k)

     ! Case 1 !
     if(s(i,j,k)*s(i+1,j,k) .le. 0.d0) then
     
       if (abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k))) .gt. tol) then

       th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k)))
       Tx_plus = (ht_Tsat-T_o(i,j,1))/th + Tij

       else

       !th = abs(s(i-1,j,k))/(abs(s(i-1,j,k))+abs(s(i+1,j,k)))
       !Tx_plus = (ht_Tsat-T_o(i-1,j,1))/th + T_o(i-1,j,1)
       !Tx_plus = (ht_Tsat-T_o(i-1,j,1))/th + Tij

       th = tol
       Tx_plus = (ht_Tsat-T_o(i,j,1))/th + Tij  

       end if 
     end if
     ! End of Case 1 !


     ! Case 2 !
     if(s(i,j,k)*s(i-1,j,k).le.0.d0) then
    
       if(abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i-1,j,k))) .gt. tol) then

       th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i-1,j,k)))
       Tx_mins = (ht_Tsat-T_o(i,j,1))/th + Tij

       else

       !th = abs(s(i+1,j,k))/(abs(s(i+1,j,k))+abs(s(i-1,j,k)))
       !Tx_mins = (ht_Tsat-T_o(i+1,j,1))/th + T_o(i+1,j,1)
       !Tx_mins = (ht_Tsat-T_o(i+1,j,1))/th + Tij
 
       th = tol
       Tx_mins = (ht_Tsat-T_o(i,j,1))/th + Tij

       
       end if
     end if
     ! End of Case 2 !

    ! Case 3 !
    if(s(i,j,k)*s(i,j+1,k).le.0.d0) then

      if (abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k))) .gt. tol) then

      th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k)))
      Ty_plus = (ht_Tsat-T_o(i,j,1))/th + Tij

      else

      !th = abs(s(i,j-1,k))/(abs(s(i,j-1,k))+abs(s(i,j+1,k)))
      !Ty_plus = (ht_Tsat-T_o(i,j-1,1))/th + T_o(i,j-1,1)  
      !Ty_plus = (ht_Tsat-T_o(i,j-1,1))/th + Tij

      th = tol
      Ty_plus = (ht_Tsat-T_o(i,j,1))/th + Tij

      end if
   end if
    ! End of Case 3 !
 
    ! Case 4 !
    if(s(i,j,k)*s(i,j-1,k).le.0.d0) then

      if (abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j-1,k))) .gt. tol) then

      th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j-1,k)))
      Ty_mins = (ht_Tsat-T_o(i,j,1))/th + Tij

      else

      !th = abs(s(i,j+1,k))/(abs(s(i,j+1,k))+abs(s(i,j-1,k)))
      !Ty_mins = (ht_Tsat-T_o(i,j+1,1))/th + T_o(i,j+1,1)
      !Ty_mins = (ht_Tsat-T_o(i,j+1,1))/th + Tij

      th = tol
      Ty_mins = (ht_Tsat-T_o(i,j,1))/th + Tij      

      end if
    end if
    ! End of Case 4 !

    !Txx = (2.*Tipj/(dxp*(dxp+dxm))) + (2.*Timj/(dxm*(dxp+dxm))) - (2.*Tij/(dxp*dxm)) 
    !Tyy = (2.*Tijp/(dyp*(dyp+dym))) + (2.*Tijm/(dym*(dyp+dym))) - (2.*Tij/(dyp*dym))   

    !! Using CELL CENTER Values of alpha

goto 100

    alphax_plus = (thco(i,j,k)/cp(i,j,k)+thco(i+1,j,k)/cp(i+1,j,k))*0.5!*(mph_thco2/mph_cp2)
    alphax_mins = (thco(i,j,k)/cp(i,j,k)+thco(i-1,j,k)/cp(i-1,j,k))*0.5!*(mph_thco2/mph_cp2)
    alphay_plus = (thco(i,j,k)/cp(i,j,k)+thco(i,j+1,k)/cp(i,j+1,k))*0.5!*(mph_thco2/mph_cp2)
    alphay_mins = (thco(i,j,k)/cp(i,j,k)+thco(i,j-1,k)/cp(i,j-1,k))*0.5!*(mph_thco2/mph_cp2)

    !! Conditions required when using CELL CENTER Values of alpha

    if (s(i,j,k)*s(i+1,j,k) .le. 0.d0) then

        th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k)))

        alphax_plus = (thco(i,j,k)/cp(i,j,k))*th + &
                          (1-th)*(thco(i+1,j,k)/cp(i+1,j,k))

    end if

    if (s(i,j,k)*s(i-1,j,k) .le. 0.d0) then

        th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i-1,j,k)))

        alphax_mins = (thco(i,j,k)/cp(i,j,k))*th + &
                          (1.-th)*(thco(i-1,j,k)/cp(i-1,j,k)) 

    end if

    if (s(i,j,k)*s(i,j+1,k) .le. 0.d0) then

        th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k)))

        alphay_plus = (thco(i,j,k)/cp(i,j,k))*th + &
                          (1.-th)*(thco(i,j+1,k)/cp(i,j+1,k))

    end if

    if (s(i,j,k)*s(i,j-1,k) .le. 0.d0) then

        th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j-1,k)))

        alphay_mins = (thco(i,j,k)/cp(i,j,k))*th + &
                          (1.-th)*(thco(i,j-1,k)/cp(i,j-1,k))

    end if  

100 continue

!_____________NON DIMENSIONAL FORM__________________________________!

    !alphax_plus = alphax_plus*(mph_thco2/mph_cp2)
    !alphax_mins = alphax_mins*(mph_thco2/mph_cp2)
    !alphay_plus = alphay_plus*(mph_thco2/mph_cp2)
    !alphay_mins = alphay_mins*(mph_thco2/mph_cp2)

!_____________NON DIMENSIONAL FORM__________________________________!

    alphax_plus = (thco(i,j,k)/cp(i,j,k))*(inRe/ht_Pr)
    alphax_mins = (thco(i,j,k)/cp(i,j,k))*(inRe/ht_Pr)
    alphay_plus = (thco(i,j,k)/cp(i,j,k))*(inRe/ht_Pr)
    alphay_mins = (thco(i,j,k)/cp(i,j,k))*(inRe/ht_Pr)

    !alphax_plus = alphax_plus*(inRe/ht_Pr)
    !alphax_mins = alphax_mins*(inRe/ht_Pr)
    !alphay_plus = alphay_plus*(inRe/ht_Pr)
    !alphay_mins = alphay_mins*(inRe/ht_Pr)
 
!_____________RHS TERM_______________________________________________!

    Txx = (alphax_plus*(Tx_plus-Tij)/dx - alphax_mins*(Tij-Tx_mins)/dx)/dx
    Tyy = (alphay_plus*(Ty_plus-Tij)/dy - alphay_mins*(Tij-Ty_mins)/dy)/dy

    T_rhs(i,j,k) = ((-(u_plus*(Tij-Tx_mins)/dx+u_mins*(Tx_plus-Tij)/dx)&
                     -(v_plus*(Tij-Ty_mins)/dy+v_mins*(Ty_plus-Tij)/dy))&
                     +(Txx +Tyy))      

     end do
  end do 

end subroutine Heat_RHS_3D
