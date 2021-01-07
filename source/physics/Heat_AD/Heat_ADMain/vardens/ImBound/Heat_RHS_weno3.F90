subroutine Heat_RHS_weno3(T_rhs, T_o, T_g, u, v, dx, dy, dz,inRe, ix1,ix2, jy1,jy2,&
                    rho1x,rho2x,rho1y,rho2y,alph,pf,s,mdot,nrmx,nrmy,smrh,curv,lambda,phip)

  use Heat_AD_data
  use Multiphase_data, only: mph_cp2,mph_thco2, mph_rho2,mph_rho1
  use Heat_AD_interface, only: Heat_GFMstencil_o1, Heat_GFMstencil_o2
  use Driver_data, only: dr_dt
 
#include "Heat_AD.h"

  implicit none
  real, dimension(:,:,:), intent(inout) :: T_rhs
  real, dimension(:,:,:), intent(in) :: T_o, T_g
  real, dimension(:,:,:), intent(in) :: u,v
  real, intent(in) :: dx, dy, dz, inRe
  integer, intent(in) :: ix1, ix2, jy1, jy2
  real, dimension(:,:,:),intent(in) :: rho1x,rho2x,rho1y,rho2y,alph
  real, dimension(:,:,:),intent(in) :: pf,s,mdot,nrmx,nrmy,smrh,curv,lambda,phip

  real :: T_res

  integer :: i,j,k

  real :: Tx_plus, Tx_mins, Ty_plus, Ty_mins, Tij
  real :: ul, ur, vl, vr
  real :: Txx, Tyy
  real :: coeff
  real :: tol
  real :: rhoc,rhoxm,rhoym,rhoxp,rhoyp
  real :: thxp1, thxm1
  real :: thyp1, thym1
  real :: thxp2, thxm2
  real :: thyp2, thym2

  real :: eps, &
          s1r,s2r,s3r,s4r,s5r,s1l,s2l,s3l,s4l,s5l, &
          rIS1r,rIS2r,rIS3r,rIS1l,rIS2l,rIS3l, &
          aT1r,aT2r,aT3r,aT1l,aT2l,aT3l, &
          a1r,a2r,a3r,a1l,a2l,a3l, &
          fT1r,fT2r,fT3r,fT1l,fT2l,fT3l, &
          frx,flx,fry,fly

  tol = 0.01
  eps = 1E-15

  k = 1

  do j=jy1,jy2
     do i=ix1,ix2

     ul = u(i,j,k)
     ur = u(i+1,j,k)
     vl = v(i,j,k)
     vr = v(i,j+1,k)

     coeff = inRe/ht_Pr

     Tx_plus = T_o(i+1,j,k)
     Tx_mins = T_o(i-1,j,k)

     Ty_plus = T_o(i,j+1,k)
     Ty_mins = T_o(i,j-1,k)

     Tij = T_o(i,j,k)

     thxp1 = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k)))
     thxm1 = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i-1,j,k)))
     thyp1 = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k)))
     thym1 = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j-1,k)))

     !______________________Diffusion Terms_______________________!

     ! Case 1 !
     if(s(i,j,k)*s(i+1,j,k) .le. 0.d0 .and. lambda(i,j,k) .lt. 0.0 .and. lambda(i+1,j,k) .lt. 0.0)&
     call Heat_GFMstencil_o1(Tx_plus,Tij,ht_Tsat,max(tol,thxp1))
     ! End of Case 1 !

     ! Case 2 !
     if(s(i,j,k)*s(i-1,j,k) .le. 0.d0 .and. lambda(i,j,k) .lt. 0.0 .and. lambda(i-1,j,k) .lt. 0.0) &
     call Heat_GFMstencil_o1(Tx_mins,Tij,ht_Tsat,max(tol,thxm1))
     ! End of Case 2 !

     ! Case 3 !
     if(s(i,j,k)*s(i,j+1,k) .le. 0.d0 .and. lambda(i,j,k) .lt. 0.0 .and. lambda(i,j+1,k) .lt. 0.0) &
     call Heat_GFMstencil_o1(Ty_plus,Tij,ht_Tsat,max(tol,thyp1))
     ! End of Case 3 !

     ! Case 4 !
     if(s(i,j,k)*s(i,j-1,k) .le. 0.d0 .and. lambda(i,j,k) .lt. 0.0 .and. lambda(i,j-1,k) .lt. 0.0) &
     call Heat_GFMstencil_o1(Ty_mins,Tij,ht_Tsat,max(tol,thym1))
     ! End of Case 4 ! 
 
     !______________________IB Terms_______________________!

     ! Case 1 !
     if(lambda(i,j,k)*lambda(i+1,j,k) .le. 0.d0) Tx_plus = T_g(i+1,j,k)
     ! End of Case 1 !

     ! Case 2 !
     if(lambda(i,j,k)*lambda(i-1,j,k) .le. 0.d0) Tx_mins = T_g(i-1,j,k)
     ! End of Case 2 !

     ! Case 3 !
     if(lambda(i,j,k)*lambda(i,j+1,k) .le. 0.d0) Ty_plus = T_g(i,j+1,k)
     ! End of Case 3 !

     ! Case 4 !
     if(lambda(i,j,k)*lambda(i,j-1,k) .le. 0.d0) Ty_mins = T_g(i,j-1,k)
     ! End of Case 4 !

     !______________________Advection Terms_______________________!

     !----------------- WENO3 X-Direction ------------!
     if (ur .gt. 0) then     ! u = (+) Downwind

        s1r = T_o(i-2,j,k)
        s2r = T_o(i-1,j,k)
        s3r = T_o(i,j,k)
        s4r = T_o(i+1,j,k)
        s5r = T_o(i+2,j,k)

        rIS1r = 13./12.*(    s1r  - 2.*s2r +    s3r )**2. &
              +  1./4. *(    s1r  - 4.*s2r + 3.*s3r )**2.
        rIS2r = 13./12.*(    s2r  - 2.*s3r +    s4r )**2. &
              +  1./4. *(    s2r           -    s4r )**2.
        rIS3r = 13./12.*(    s3r  - 2.*s4r +    s5r )**2. &
              +  1./4. *( 3.*s3r  - 4.*s4r +    s5r )**2.

        aT1r = 1./10. /  ( eps + rIS1r )**2.
        aT2r = 6./10. /  ( eps + rIS2r )**2.
        aT3r = 3./10. /  ( eps + rIS3r )**2.

        a1r = aT1r / ( aT1r + aT2r +aT3r )
        a2r = aT2r / ( aT1r + aT2r +aT3r )
        a3r = aT3r / ( aT1r + aT2r +aT3r )

        fT1r =  2./6.*s1r - 7./6.*s2r + 11./6.*s3r
        fT2r = -1./6.*s2r + 5./6.*s3r +  2./6.*s4r
        fT3r =  2./6.*s3r + 5./6.*s4r -  1./6.*s5r

      else                  ! u = (-) Upwind

        s1r = T_o(i-1,j,k)
        s2r = T_o(i,j,k)
        s3r = T_o(i+1,j,k)
        s4r = T_o(i+2,j,k)
        s5r = T_o(i+3,j,k)

        rIS1r = 13./12.*(    s1r  - 2.*s2r +    s3r )**2. &
              +  1./4. *(    s1r  - 4.*s2r + 3.*s3r )**2.
        rIS2r = 13./12.*(    s2r  - 2.*s3r +    s4r )**2. &
              +  1./4. *(    s2r           -    s4r )**2.
        rIS3r = 13./12.*(    s3r  - 2.*s4r +    s5r )**2. &
              +  1./4. *( 3.*s3r  - 4.*s4r +    s5r )**2.

        aT1r = 3./10. /  ( eps + rIS1r )**2.
        aT2r = 6./10. /  ( eps + rIS2r )**2.
        aT3r = 1./10. /  ( eps + rIS3r )**2.

        a1r = aT1r / ( aT1r + aT2r +aT3r )
        a2r = aT2r / ( aT1r + aT2r +aT3r )
        a3r = aT3r / ( aT1r + aT2r +aT3r )

        fT1r = -1./6.*s1r + 5./6.*s2r +  2./6.*s3r
        fT2r =  2./6.*s2r + 5./6.*s3r -  1./6.*s4r
        fT3r =  11./6.*s3r - 7./6.*s4r + 2./6.*s5r

     end if

     if (ul .gt. 0) then     ! u = (+) Downwind  

        s1l = T_o(i-3,j,k)
        s2l = T_o(i-2,j,k)
        s3l = T_o(i-1,j,k)
        s4l = T_o(i,j,k)
        s5l = T_o(i+1,j,k)

        rIS1l = 13./12.*(    s1l  - 2.*s2l +    s3l )**2. &
              +  1./4. *(    s1l  - 4.*s2l + 3.*s3l )**2.
        rIS2l = 13./12.*(    s2l  - 2.*s3l +    s4l )**2. &
              +  1./4. *(    s2l           -    s4l )**2.
        rIS3l = 13./12.*(    s3l  - 2.*s4l +    s5l )**2. &
              +  1./4. *( 3.*s3l  - 4.*s4l +    s5l )**2.

        aT1l = 1./10. /  ( eps + rIS1l )**2.
        aT2l = 6./10. /  ( eps + rIS2l )**2.
        aT3l = 3./10. /  ( eps + rIS3l )**2.

        a1l = aT1l / ( aT1l + aT2l +aT3l )
        a2l = aT2l / ( aT1l + aT2l +aT3l )
        a3l = aT3l / ( aT1l + aT2l +aT3l )

        fT1l =  2./6.*s1l - 7./6.*s2l + 11./6.*s3l
        fT2l = -1./6.*s2l + 5./6.*s3l +  2./6.*s4l
        fT3l =  2./6.*s3l + 5./6.*s4l -  1./6.*s5l

     else                   ! u = (-) Upwind

        s1l = T_o(i-2,j,k)
        s2l = T_o(i-1,j,k)
        s3l = T_o(i,j,k)
        s4l = T_o(i+1,j,k)
        s5l = T_o(i+2,j,k)

        rIS1l = 13./12.*(    s1l  - 2.*s2l +    s3l )**2. &
              +  1./4. *(    s1l  - 4.*s2l + 3.*s3l )**2.
        rIS2l = 13./12.*(    s2l  - 2.*s3l +    s4l )**2. &
              +  1./4. *(    s2l           -    s4l )**2.
        rIS3l = 13./12.*(    s3l  - 2.*s4l +    s5l )**2. &
              +  1./4. *( 3.*s3l  - 4.*s4l +    s5l )**2.

        aT1l = 3./10. /  ( eps + rIS1l )**2.
        aT2l = 6./10. /  ( eps + rIS2l )**2.
        aT3l = 1./10. /  ( eps + rIS3l )**2.

        a1l = aT1l / ( aT1l + aT2l +aT3l )
        a2l = aT2l / ( aT1l + aT2l +aT3l )
        a3l = aT3l / ( aT1l + aT2l +aT3l )

        fT1l = -1./6.*s1l + 5./6.*s2l +  2./6.*s3l
        fT2l =  2./6.*s2l + 5./6.*s3l -  1./6.*s4l
        fT3l =  11./6.*s3l - 7./6.*s4l + 2./6.*s5l

     end if

     !---------------------------------------------------------
     !- WENO3 interpolated TEMP values at cell face
     !---------------------------------------------------------
     frx = a1r*fT1r + a2r*fT2r + a3r*fT3r
     flx = a1l*fT1l + a2l*fT2l + a3l*fT3l
     !---------------------------------------------------------
     !---------------------------------------------------------

     !----------------- WENO3 Y-Direction ------------!
     if (vr .gt. 0) then     ! u = (+) Downwind

        s1r = T_o(i,j-2,k)
        s2r = T_o(i,j-1,k)
        s3r = T_o(i,j,k)
        s4r = T_o(i,j+1,k)
        s5r = T_o(i,j+2,k)

        rIS1r = 13./12.*(    s1r  - 2.*s2r +    s3r )**2. &
              +  1./4. *(    s1r  - 4.*s2r + 3.*s3r )**2.
        rIS2r = 13./12.*(    s2r  - 2.*s3r +    s4r )**2. &
              +  1./4. *(    s2r           -    s4r )**2.
        rIS3r = 13./12.*(    s3r  - 2.*s4r +    s5r )**2. &
              +  1./4. *( 3.*s3r  - 4.*s4r +    s5r )**2.

        aT1r = 1./10. /  ( eps + rIS1r )**2.
        aT2r = 6./10. /  ( eps + rIS2r )**2.
        aT3r = 3./10. /  ( eps + rIS3r )**2.

        a1r = aT1r / ( aT1r + aT2r +aT3r )
        a2r = aT2r / ( aT1r + aT2r +aT3r )
        a3r = aT3r / ( aT1r + aT2r +aT3r )

        fT1r =  2./6.*s1r - 7./6.*s2r + 11./6.*s3r
        fT2r = -1./6.*s2r + 5./6.*s3r +  2./6.*s4r
        fT3r =  2./6.*s3r + 5./6.*s4r -  1./6.*s5r

     else                   ! u = (-) Upwind

        s1r = T_o(i,j-1,k)
        s2r = T_o(i,j,k)
        s3r = T_o(i,j+1,k)
        s4r = T_o(i,j+2,k)
        s5r = T_o(i,j+3,k)

        rIS1r = 13./12.*(    s1r  - 2.*s2r +    s3r )**2. &
              +  1./4. *(    s1r  - 4.*s2r + 3.*s3r )**2.
        rIS2r = 13./12.*(    s2r  - 2.*s3r +    s4r )**2. &
              +  1./4. *(    s2r           -    s4r )**2.
        rIS3r = 13./12.*(    s3r  - 2.*s4r +    s5r )**2. &
              +  1./4. *( 3.*s3r  - 4.*s4r +    s5r )**2.

        aT1r = 3./10. /  ( eps + rIS1r )**2.
        aT2r = 6./10. /  ( eps + rIS2r )**2.
        aT3r = 1./10. /  ( eps + rIS3r )**2.

        a1r = aT1r / ( aT1r + aT2r +aT3r )
        a2r = aT2r / ( aT1r + aT2r +aT3r )
        a3r = aT3r / ( aT1r + aT2r +aT3r )

        fT1r = -1./6.*s1r + 5./6.*s2r +  2./6.*s3r
        fT2r =  2./6.*s2r + 5./6.*s3r -  1./6.*s4r
        fT3r =  11./6.*s3r - 7./6.*s4r + 2./6.*s5r

     end if

     if (vl .gt. 0) then     ! u = (+) Downwind

        s1l = T_o(i,j-3,k)
        s2l = T_o(i,j-2,k)
        s3l = T_o(i,j-1,k)
        s4l = T_o(i,j,k)
        s5l = T_o(i,j+1,k)

        rIS1l = 13./12.*(    s1l  - 2.*s2l +    s3l )**2. &
              +  1./4. *(    s1l  - 4.*s2l + 3.*s3l )**2.
        rIS2l = 13./12.*(    s2l  - 2.*s3l +    s4l )**2. &
              +  1./4. *(    s2l           -    s4l )**2.
        rIS3l = 13./12.*(    s3l  - 2.*s4l +    s5l )**2. &
              +  1./4. *( 3.*s3l  - 4.*s4l +    s5l )**2.

        aT1l = 1./10. /  ( eps + rIS1l )**2.
        aT2l = 6./10. /  ( eps + rIS2l )**2.
        aT3l = 3./10. /  ( eps + rIS3l )**2.

        a1l = aT1l / ( aT1l + aT2l +aT3l )
        a2l = aT2l / ( aT1l + aT2l +aT3l )
        a3l = aT3l / ( aT1l + aT2l +aT3l )

        fT1l =  2./6.*s1l - 7./6.*s2l + 11./6.*s3l
        fT2l = -1./6.*s2l + 5./6.*s3l +  2./6.*s4l
        fT3l =  2./6.*s3l + 5./6.*s4l -  1./6.*s5l

     else                    ! u = (-) Upwind

        s1l = T_o(i,j-2,k)
        s2l = T_o(i,j-1,k)
        s3l = T_o(i,j,k)
        s4l = T_o(i,j+1,k)
        s5l = T_o(i,j+2,k)

        rIS1l = 13./12.*(    s1l  - 2.*s2l +    s3l )**2. &
              +  1./4. *(    s1l  - 4.*s2l + 3.*s3l )**2.
        rIS2l = 13./12.*(    s2l  - 2.*s3l +    s4l )**2. &
              +  1./4. *(    s2l           -    s4l )**2.
        rIS3l = 13./12.*(    s3l  - 2.*s4l +    s5l )**2. &
              +  1./4. *( 3.*s3l  - 4.*s4l +    s5l )**2.

        aT1l = 3./10. /  ( eps + rIS1l )**2.
        aT2l = 6./10. /  ( eps + rIS2l )**2.
        aT3l = 1./10. /  ( eps + rIS3l )**2.

        a1l = aT1l / ( aT1l + aT2l +aT3l )
        a2l = aT2l / ( aT1l + aT2l +aT3l )
        a3l = aT3l / ( aT1l + aT2l +aT3l )

        fT1l = -1./6.*s1l + 5./6.*s2l +  2./6.*s3l
        fT2l =  2./6.*s2l + 5./6.*s3l -  1./6.*s4l
        fT3l =  11./6.*s3l - 7./6.*s4l + 2./6.*s5l

     end if

     !---------------------------------------------------------
     !- WENO3 interpolated TEMP values at cell face
     !---------------------------------------------------------
     fry = a1r*fT1r + a2r*fT2r + a3r*fT3r
     fly = a1l*fT1l + a2l*fT2l + a3l*fT3l
     !---------------------------------------------------------
     !---------------------------------------------------------

!_______________________________RHS TERM______________________________________!

    Txx = alph(i,j,k)*(coeff*(Tx_plus-Tij)/dx - coeff*(Tij-Tx_mins)/dx)/dx
    Tyy = alph(i,j,k)*(coeff*(Ty_plus-Tij)/dy - coeff*(Tij-Ty_mins)/dy)/dy

    if(lambda(i,j,k) .lt. 0.0) then

    T_rhs(i,j,k) = - (frx*ur - flx*ul)/dx &
                   - (fry*vr - fly*vl)/dy &
                   + Txx + Tyy

    else

    T_rhs(i,j,k) = (1.0 - Tij)/dr_dt !Txx + Tyy

    end if

    end do
  end do 

end subroutine Heat_RHS_weno3
