subroutine Heat_extrapGradT_weno3(Tnl,Tnv,T,s,pf,dx,dy,dz,nx,ny,ix1,ix2,jy1,jy2,Tnl_res,Tnv_res,mflg)

#include "Flash.h"

   implicit none
   real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
   real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny,mflg
   real, intent(in) :: dx,dy,dz
   integer, intent(in) :: ix1,ix2,jy1,jy2
   real, intent(out) :: Tnl_res,Tnv_res

   integer :: i,j,k

   real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: Tnl_o,Tnv_o,Tnl_i,Tnv_i

   real :: dt_ext,maxN

   integer :: l_counter, v_counter, step
   logical :: int_xp, int_xm, int_yp, int_ym

   real :: nxl,nxr,nyl,nyr

   real :: eps, &
           s1r,s2r,s3r,s4r,s5r,s1l,s2l,s3l,s4l,s5l, &
           rIS1r,rIS2r,rIS3r,rIS1l,rIS2l,rIS3l, &
           aT1r,aT2r,aT3r,aT1l,aT2l,aT3l, &
           a1r,a2r,a3r,a1l,a2l,a3l, &
           fT1r,fT2r,fT3r,fT1l,fT2l,fT3l, &
           frx,flx,fry,fly

   l_counter = 0
   v_counter = 0
   k = 1

   eps = 1E-15

   maxN = max(maxval(abs(nx)),maxval(abs(ny)))

   dt_ext = 0.05d0*(dx/maxN)

   Tnl_i = Tnl
   Tnv_i = Tnv

   Tnl_o = Tnl

   Tnv_o = Tnv

   do j=jy1,jy2
     do i=ix1,ix2

     nxl = sign(1.,s(i,j,k))*(nx(i,j,k)+nx(i-1,j,k))*0.5
     nxr = sign(1.,s(i,j,k))*(nx(i,j,k)+nx(i+1,j,k))*0.5

     nyl = sign(1.,s(i,j,k))*(ny(i,j,k)+ny(i,j-1,k))*0.5
     nyr = sign(1.,s(i,j,k))*(ny(i,j,k)+ny(i,j+1,k))*0.5

     !----------------- WENO3 X-Direction ------------!
     if (nxr .gt. 0) then     ! nx = (+) Downwind

        if(pf(i,j,k) .eq. 0) then
                s1r = Tnv_o(i-2,j,k)
                s2r = Tnv_o(i-1,j,k)
                s3r = Tnv_o(i,j,k)
                s4r = Tnv_o(i+1,j,k)
                s5r = Tnv_o(i+2,j,k)
        else
                s1r = Tnl_o(i-2,j,k)
                s2r = Tnl_o(i-1,j,k)
                s3r = Tnl_o(i,j,k)
                s4r = Tnl_o(i+1,j,k)
                s5r = Tnl_o(i+2,j,k)
        end if

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

      else                  ! nx = (-) Upwind

        if(pf(i,j,k) .eq. 0) then
                s1r = Tnv_o(i-1,j,k)
                s2r = Tnv_o(i,j,k)
                s3r = Tnv_o(i+1,j,k)
                s4r = Tnv_o(i+2,j,k)
                s5r = Tnv_o(i+3,j,k)
        else
                s1r = Tnl_o(i-1,j,k)
                s2r = Tnl_o(i,j,k)
                s3r = Tnl_o(i+1,j,k)
                s4r = Tnl_o(i+2,j,k)
                s5r = Tnl_o(i+3,j,k)
        end if

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
        fT3r = 11./6.*s3r - 7./6.*s4r + 2./6.*s5r

     end if

     if (nxl .gt. 0) then     ! nx = (+) Downwind  

        if(pf(i,j,k) .eq. 0) then
                s1l = Tnv_o(i-3,j,k)
                s2l = Tnv_o(i-2,j,k)
                s3l = Tnv_o(i-1,j,k)
                s4l = Tnv_o(i,j,k)
                s5l = Tnv_o(i+1,j,k)
        else
                s1l = Tnl_o(i-3,j,k)
                s2l = Tnl_o(i-2,j,k)
                s3l = Tnl_o(i-1,j,k)
                s4l = Tnl_o(i,j,k)
                s5l = Tnl_o(i+1,j,k)
        end if

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

     else                   ! nx = (-) Upwind

        if(pf(i,j,k) .eq. 0) then
                s1l = Tnv_o(i-2,j,k)
                s2l = Tnv_o(i-1,j,k)
                s3l = Tnv_o(i,j,k)
                s4l = Tnv_o(i+1,j,k)
                s5l = Tnv_o(i+2,j,k)
        else
                s1l = Tnl_o(i-2,j,k)
                s2l = Tnl_o(i-1,j,k)
                s3l = Tnl_o(i,j,k)
                s4l = Tnl_o(i+1,j,k)
                s5l = Tnl_o(i+2,j,k)
        end if

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
        fT3l = 11./6.*s3l - 7./6.*s4l + 2./6.*s5l

     end if

     !---------------------------------------------------------
     !- WENO3 interpolated HEAT FLUX values at cell face
     !---------------------------------------------------------
     frx = a1r*fT1r + a2r*fT2r + a3r*fT3r
     flx = a1l*fT1l + a2l*fT2l + a3l*fT3l
     !---------------------------------------------------------

     !----------------- WENO3 Y-Direction ------------!
     if (nyr .gt. 0) then     ! ny = (+) Downwind

        if(pf(i,j,k) .eq. 0) then
                s1r = Tnv_o(i,j-2,k)
                s2r = Tnv_o(i,j-1,k)
                s3r = Tnv_o(i,j,k)
                s4r = Tnv_o(i,j+1,k)
                s5r = Tnv_o(i,j+2,k)
        else
                s1r = Tnl_o(i,j-2,k)
                s2r = Tnl_o(i,j-1,k)
                s3r = Tnl_o(i,j,k)
                s4r = Tnl_o(i,j+1,k)
                s5r = Tnl_o(i,j+2,k)
        end if

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

     else                   ! ny = (-) Upwind

        if(pf(i,j,k) .eq. 0) then
                s1r = Tnv_o(i,j-1,k)
                s2r = Tnv_o(i,j,k)
                s3r = Tnv_o(i,j+1,k)
                s4r = Tnv_o(i,j+2,k)
                s5r = Tnv_o(i,j+3,k)
        else  
                s1r = Tnl_o(i,j-1,k)
                s2r = Tnl_o(i,j,k)
                s3r = Tnl_o(i,j+1,k)
                s4r = Tnl_o(i,j+2,k)
                s5r = Tnl_o(i,j+3,k)
        end if

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
        fT3r = 11./6.*s3r - 7./6.*s4r +  2./6.*s5r

     end if

     if (nyl .gt. 0) then     ! ny = (+) Downwind

        if(pf(i,j,k) .eq. 0) then
                s1l = Tnv_o(i,j-3,k)
                s2l = Tnv_o(i,j-2,k)
                s3l = Tnv_o(i,j-1,k)
                s4l = Tnv_o(i,j,k)
                s5l = Tnv_o(i,j+1,k)
        else  
                s1l = Tnl_o(i,j-3,k)
                s2l = Tnl_o(i,j-2,k)
                s3l = Tnl_o(i,j-1,k)
                s4l = Tnl_o(i,j,k)
                s5l = Tnl_o(i,j+1,k)
        end if

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

     else                    ! ny = (-) Upwind

        if(pf(i,j,k) .eq. 0) then
                s1l = Tnv_o(i,j-2,k)
                s2l = Tnv_o(i,j-1,k)
                s3l = Tnv_o(i,j,k)
                s4l = Tnv_o(i,j+1,k)
                s5l = Tnv_o(i,j+2,k)
        else  
                s1l = Tnl_o(i,j-2,k)
                s2l = Tnl_o(i,j-1,k)
                s3l = Tnl_o(i,j,k)
                s4l = Tnl_o(i,j+1,k)
                s5l = Tnl_o(i,j+2,k)
        end if

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
     !- WENO3 interpolated HEAT FLUX values at cell face
     !---------------------------------------------------------
     fry = a1r*fT1r + a2r*fT2r + a3r*fT3r
     fly = a1l*fT1l + a2l*fT2l + a3l*fT3l
     !---------------------------------------------------------

     Tnl(i,j,k) = Tnl_i(i,j,k) + dt_ext*pf(i,j,k)*(-(frx*nxr - flx*nxl)/dx &
                                                   -(fry*nyr - fly*nyl)/dy )

     Tnv(i,j,k) = Tnv_i(i,j,k) + dt_ext*(1.0-pf(i,j,k))*(-(frx*nxr - flx*nxl)/dx &
                                                         -(fry*nyr - fly*nyl)/dy )
     end do
   end do   

   Tnl_res = sum(sum((mflg(ix1:ix2,jy1:jy2,k)*(Tnl_o(ix1:ix2,jy1:jy2,k)-Tnl(ix1:ix2,jy1:jy2,k)))**2,1),1)
   Tnv_res = sum(sum((mflg(ix1:ix2,jy1:jy2,k)*(Tnv_o(ix1:ix2,jy1:jy2,k)-Tnv(ix1:ix2,jy1:jy2,k)))**2,1),1)

   Tnl_res = sqrt(Tnl_res/size(Tnl(ix1:ix2,jy1:jy2,k)))
   Tnv_res = sqrt(Tnv_res/size(Tnv(ix1:ix2,jy1:jy2,k)))

end subroutine Heat_extrapGradT_weno3
