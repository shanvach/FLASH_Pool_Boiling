SUBROUTINE ins_rhs2d_weno3(uni,vni,ru1,ix1,ix2,jy1,jy2,dx,dy,ru,rv, &
                           visc,rho1x,rho2x,rho1y,rho2y,gravX,gravY, &
                           mdot,smrh,xnorm,ynorm,crv,s,pf,temp,blockID)

  !***************************************************************
  ! This subroutine computes the discretization of the RHS of the 
  ! Helmholtz equation on a staggered uniform grid.
  !
  ! Input:  uni,vni     = velocity at timestep n
  !         ru1         = molecular viscosity !- kpd - Inverse Reynolds No
  !         ix1,ix2     = starting and ending x indices
  !         jy1,jy2     = starting and ending y indices
  !         dx,dy       = grid spacing in x and y directions
  !
  ! Output: ru,rv    = u and v momentum for Helmholtz RHS
  !**************************************************************

      use Driver_interface, ONLY : Driver_abortFlash

      use Driver_data,      ONLY : dr_nstep

      use RuntimeParameters_interface, ONLY : RuntimeParameters_get

      use IncompNS_data, ONLY : ins_iConvU, ins_convvel
 
      use Multiphase_data, ONLY: mph_rho2,mph_sten

      use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_getBlkCenterCoords, Grid_getDeltas

      use Simulation_data, only: sim_yMax, sim_yMin

#include "Flash.h"
#include "constants.h"

      implicit none
      INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2
      REAL, INTENT(IN):: ru1, dx, dy, gravX,gravY
      REAL, DIMENSION(:,:,:), INTENT(IN):: uni, vni, visc, rho1x, rho2x, rho1y, rho2y
      REAL, DIMENSION(:,:,:), INTENT(IN) :: xnorm,ynorm,mdot,smrh,crv,s,pf,temp
      REAL, DIMENSION(:,:,:), INTENT(OUT):: ru, rv
      INTEGER, INTENT(IN) :: blockID

      INTEGER:: i, j
      REAL:: dx1, dy1, Mdens, th, cri, crc
      ! x-component variables
      REAL:: uxplus, uxminus, vxplus, vxminus, wxplus, wxminus, &
             uyplus, uyminus
      REAL:: dudxp, dudxm, dudyp, dudym, dvdxp, dvdxm
      REAL:: dsdxp, dsdxm, dsdyp, dsdym
      REAL:: tvjp, tvjm
      REAL:: txxp, txxm, tyyp, tyym
      REAL:: txyp, txym
      ! new y-component variables
      REAL:: vyplus, vyminus
      REAL:: dvdyp, dvdym
      REAL:: tvip, tvim
      INTEGER, parameter :: kz1 = 1

      real :: eps, &
              s1r,s2r,s3r,s4r,s5r,s1l,s2l,s3l,s4l,s5l, &
              rIS1r,rIS2r,rIS3r,rIS1l,rIS2l,rIS3l, &
              aT1r,aT2r,aT3r,aT1l,aT2l,aT3l, &
              a1r,a2r,a3r,a1l,a2l,a3l, &
              fT1r,fT2r,fT3r,fT1l,fT2l,fT3l, &
              frx,flx,fry,fly

      real :: aicc, aixl, aixr, aiyl, aiyr

      real :: ycell,Lb
      real :: del(MDIM),bsize(MDIM),coord(MDIM)
      real, dimension(2,MDIM) :: boundBox

      ! grid spacings
      dx1 = 1.0/dx
      dy1 = 1.0/dy

      eps = 1E-15
  
      Lb  = sim_yMax-sim_yMin-2.0

      call Grid_getDeltas(blockID,del)
      call Grid_getBlkCenterCoords(blockId,coord)
      call Grid_getBlkBoundBox(blockId,boundBox)

      bsize(:) = boundBox(2,:) - boundBox(1,:)

      !++++++++++  U-COMPONENT  ++++++++++
       do j = jy1,jy2
          do i = ix1,ix2+1

             ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                      real(j - NGUARD - 1)*del(JAXIS)  +  &
                      0.5*del(JAXIS)
  
             ! Advection Terms

             uxplus = (uni(i+1,j,kz1) + uni(i,j,kz1))*0.5 
             uxminus = (uni(i-1,j,kz1) + uni(i,j,kz1))*0.5

             vyplus = (vni(i,j+1,kz1) + vni(i-1,j+1,kz1))*0.5
             vyminus = (vni(i,j,kz1) + vni(i-1,j,kz1))*0.5

             uyplus = (uni(i,j+1,kz1) + uni(i,j,kz1))*0.5
             uyminus = (uni(i,j,kz1) + uni(i,j-1,kz1))*0.5

             !----------------- WENO3 X-Direction ------------!
             if (uxplus .gt. 0) then     ! u = (+) Downwind

                s1r = uni(i-2,j,kz1)
                s2r = uni(i-1,j,kz1)
                s3r = uni(i,j,kz1)
                s4r = uni(i+1,j,kz1)
                s5r = uni(i+2,j,kz1)

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

                s1r = uni(i-1,j,kz1)
                s2r = uni(i,j,kz1)
                s3r = uni(i+1,j,kz1)
                s4r = uni(i+2,j,kz1)
                s5r = uni(i+3,j,kz1)

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


            if (uxminus .gt. 0) then     ! u = (+) Downwind  

                s1l = uni(i-3,j,kz1)
                s2l = uni(i-2,j,kz1)
                s3l = uni(i-1,j,kz1)
                s4l = uni(i,j,kz1)
                s5l = uni(i+1,j,kz1)

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

                s1l = uni(i-2,j,kz1)
                s2l = uni(i-1,j,kz1)
                s3l = uni(i,j,kz1)
                s4l = uni(i+1,j,kz1)
                s5l = uni(i+2,j,kz1)

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
          !- WENO3 interpolated VEL values at cell center
          !---------------------------------------------------------
          frx = a1r*fT1r + a2r*fT2r + a3r*fT3r
          flx = a1l*fT1l + a2l*fT2l + a3l*fT3l
          !---------------------------------------------------------
          !---------------------------------------------------------


          !----------------- WENO3 Y-Direction ------------!
          if (vyplus .gt. 0) then     ! u = (+) Downwind

                s1r = uni(i,j-2,kz1)
                s2r = uni(i,j-1,kz1)
                s3r = uni(i,j,kz1)
                s4r = uni(i,j+1,kz1)
                s5r = uni(i,j+2,kz1)

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

                s1r = uni(i,j-1,kz1)
                s2r = uni(i,j,kz1)
                s3r = uni(i,j+1,kz1)
                s4r = uni(i,j+2,kz1)
                s5r = uni(i,j+3,kz1)

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

           if (vyminus .gt. 0) then     ! u = (+) Downwind

                s1l = uni(i,j-3,kz1)
                s2l = uni(i,j-2,kz1)
                s3l = uni(i,j-1,kz1)
                s4l = uni(i,j,kz1)
                s5l = uni(i,j+1,kz1)

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

                s1l = uni(i,j-2,kz1)
                s2l = uni(i,j-1,kz1)
                s3l = uni(i,j,kz1)
                s4l = uni(i,j+1,kz1)
                s5l = uni(i,j+2,kz1)

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
             !- WENO3 interpolated VEL values at cell center
             !---------------------------------------------------------
             fry = a1r*fT1r + a2r*fT2r + a3r*fT3r
             fly = a1l*fT1l + a2l*fT2l + a3l*fT3l
             !---------------------------------------------------------
             !---------------------------------------------------------

             ! Diffusion Terms

             ! get derivatives at 1/2 locations

             ! velocity
             dudxp = (uni(i+1,j,kz1) - uni(i,j,kz1))*dx1
             dudxm = (uni(i,j,kz1)   - uni(i-1,j,kz1))*dx1
             dudyp = (uni(i,j+1,kz1) - uni(i,j,kz1))*dy1
             dudym = (uni(i,j,kz1)   - uni(i,j-1,kz1))*dy1

             ! velocity jump source terms

             ! Method 1 - Smeared Density at Cell Centers

             aicc = 0.5*(smrh(i,j,kz1)+smrh(i-1,j,kz1))*&
                    0.5*(mdot(i,j,kz1)+mdot(i-1,j,kz1))*&
                    0.5*(xnorm(i,j,kz1)+xnorm(i-1,j,kz1))

             aixr = 0.5*(smrh(i,j,kz1)+smrh(i+1,j,kz1))*&
                    0.5*(mdot(i,j,kz1)+mdot(i-1,j,kz1))*&
                    0.5*(xnorm(i,j,kz1)+xnorm(i-1,j,kz1))

             aixl = 0.5*(smrh(i-1,j,kz1)+smrh(i-2,j,kz1))*&
                    0.5*(mdot(i,j,kz1)+mdot(i-1,j,kz1))*&
                    0.5*(xnorm(i,j,kz1)+xnorm(i-1,j,kz1))

             aiyr = 0.5*(smrh(i,j+1,kz1)+smrh(i-1,j+1,kz1))*&
                    0.5*(mdot(i,j,kz1)+mdot(i-1,j,kz1))*&
                    0.5*(xnorm(i,j,kz1)+xnorm(i-1,j,kz1)) 

             aiyl = 0.5*(smrh(i,j-1,kz1)+smrh(i-1,j-1,kz1))*&
                    0.5*(mdot(i,j,kz1)+mdot(i-1,j,kz1))*&
                    0.5*(xnorm(i,j,kz1)+xnorm(i-1,j,kz1))

             ! Method 2 - Smeared Density at Cell Faces

             !aicc = (rho1x(i,j,kz1)+rho2x(i,j,kz1))*&
             !       0.5*(mdot(i,j,kz1)+mdot(i-1,j,kz1))*&
             !       0.5*(xnorm(i,j,kz1)+xnorm(i-1,j,kz1))

             !aixr = (rho1x(i+1,j,kz1)+rho2x(i+1,j,kz1))*&
             !       0.5*(mdot(i,j,kz1)+mdot(i-1,j,kz1))*&
             !       0.5*(xnorm(i,j,kz1)+xnorm(i-1,j,kz1))

             !aixl = (rho1x(i-1,j,kz1)+rho2x(i-1,j,kz1))*&
             !       0.5*(mdot(i,j,kz1)+mdot(i-1,j,kz1))*&
             !       0.5*(xnorm(i,j,kz1)+xnorm(i-1,j,kz1))

             !aiyr = (rho1x(i,j+1,kz1)+rho2x(i,j+1,kz1))*&
             !       0.5*(mdot(i,j,kz1)+mdot(i-1,j,kz1))*&
             !       0.5*(xnorm(i,j,kz1)+xnorm(i-1,j,kz1)) 

             !aiyl = (rho1x(i,j-1,kz1)+rho2x(i,j-1,kz1))*&
             !       0.5*(mdot(i,j,kz1)+mdot(i-1,j,kz1))*&
             !       0.5*(xnorm(i,j,kz1)+xnorm(i-1,j,kz1))

             dsdxp = (aixr - aicc)*dx1
             dsdxm = (aicc - aixl)*dx1
             dsdyp = (aiyr - aicc)*dy1
             dsdym = (aicc - aiyl)*dy1

             !-Variable Viscosity Implementation (ru1 is 1/Re)
             txxp = ru1*visc(i  ,j,kz1)*(dudxp+dsdxp)
             txxm = ru1*visc(i-1,j,kz1)*(dudxm+dsdxm)
             tyyp = ru1*0.25*(visc(i,j,kz1)+visc(i-1,j,kz1)+visc(i-1,j+1,kz1)+visc(i,j+1,kz1))*(dudyp+dsdyp)
             tyym = ru1*0.25*(visc(i,j,kz1)+visc(i-1,j,kz1)+visc(i-1,j-1,kz1)+visc(i,j-1,kz1))*(dudym+dsdym)

             Mdens = ( rho1x(i,j,kz1) + rho2x(i,j,kz1) )  ! Mixture inverse density

             ! calculate RHS for u-momentum
             ru(i,j,kz1) = - (frx*uxplus - flx*uxminus)/dx &
                           - (fry*vyplus - fly*vyminus)/dy &
                           + Mdens*(txxp - txxm)*dx1                 & ! diffusion - normal terms 
                           + Mdens*(tyyp - tyym)*dy1                 &
                           + gravX                                   

            !if(ycell .ge. Lb) ru(i,j,kz1) = ru(i,j,kz1) - ins_convvel(HIGH,JAXIS)*(uyplus - uyminus)*dy1 

          enddo
       enddo


    !++++++++++  V-COMPONENT  ++++++++++

       do j = jy1,jy2+1
          do i = ix1,ix2

             ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                      real(j - NGUARD - 1)*del(JAXIS)

             uxplus =  (uni(i+1,j,kz1) + uni(i+1,j-1,kz1))*0.5 
             uxminus = (uni(i,j,kz1) + uni(i,j-1,kz1))*0.5

             vyplus = (vni(i,j+1,kz1) + vni(i,j,kz1))*0.5
             vyminus = (vni(i,j,kz1) + vni(i,j-1,kz1))*0.5

             !----------------- WENO3 X-Direction ------------!
             if (uxplus .gt. 0) then     ! u = (+) Downwind

                s1r = vni(i-2,j,kz1)
                s2r = vni(i-1,j,kz1)
                s3r = vni(i,j,kz1)
                s4r = vni(i+1,j,kz1)
                s5r = vni(i+2,j,kz1)

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

                s1r = vni(i-1,j,kz1)
                s2r = vni(i,j,kz1)
                s3r = vni(i+1,j,kz1)
                s4r = vni(i+2,j,kz1)
                s5r = vni(i+3,j,kz1)

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


            if (uxminus .gt. 0) then     ! u = (+) Downwind  

                s1l = vni(i-3,j,kz1)
                s2l = vni(i-2,j,kz1)
                s3l = vni(i-1,j,kz1)
                s4l = vni(i,j,kz1)
                s5l = vni(i+1,j,kz1)

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

                s1l = vni(i-2,j,kz1)
                s2l = vni(i-1,j,kz1)
                s3l = vni(i,j,kz1)
                s4l = vni(i+1,j,kz1)
                s5l = vni(i+2,j,kz1)

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
          !- WENO3 interpolated VEL values at cell center
          !---------------------------------------------------------
          frx = a1r*fT1r + a2r*fT2r + a3r*fT3r
          flx = a1l*fT1l + a2l*fT2l + a3l*fT3l
          !---------------------------------------------------------
          !---------------------------------------------------------


          !----------------- WENO3 Y-Direction ------------!
          if (vyplus .gt. 0) then     ! u = (+) Downwind

                s1r = vni(i,j-2,kz1)
                s2r = vni(i,j-1,kz1)
                s3r = vni(i,j,kz1)
                s4r = vni(i,j+1,kz1)
                s5r = vni(i,j+2,kz1)

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

                s1r = vni(i,j-1,kz1)
                s2r = vni(i,j,kz1)
                s3r = vni(i,j+1,kz1)
                s4r = vni(i,j+2,kz1)
                s5r = vni(i,j+3,kz1)

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

           if (vyminus .gt. 0) then     ! u = (+) Downwind

                s1l = vni(i,j-3,kz1)
                s2l = vni(i,j-2,kz1)
                s3l = vni(i,j-1,kz1)
                s4l = vni(i,j,kz1)
                s5l = vni(i,j+1,kz1)

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

                s1l = vni(i,j-2,kz1)
                s2l = vni(i,j-1,kz1)
                s3l = vni(i,j,kz1)
                s4l = vni(i,j+1,kz1)
                s5l = vni(i,j+2,kz1)

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
             !- WENO3 interpolated VEL values at cell center
             !---------------------------------------------------------
             fry = a1r*fT1r + a2r*fT2r + a3r*fT3r
             fly = a1l*fT1l + a2l*fT2l + a3l*fT3l
             !---------------------------------------------------------
             !---------------------------------------------------------


             ! Diffusion Terms

             ! get derivatives at 1/2 locations
   
             ! velocity
             dvdxp = (vni(i+1,j,kz1) - vni(i,j,kz1))*dx1
             dvdxm = (vni(i,j,kz1) - vni(i-1,j,kz1))*dx1
             dvdyp = (vni(i,j+1,kz1) - vni(i,j,kz1))*dy1
             dvdym = (vni(i,j,kz1) - vni(i,j-1,kz1))*dy1

             ! velocity jump source terms

             ! Method 1 - Smeared Density at Cell Centers

             aicc = 0.5*(smrh(i,j,kz1)+smrh(i,j-1,kz1))*&
                    0.5*(mdot(i,j,kz1)+mdot(i,j-1,kz1))*&
                    0.5*(ynorm(i,j,kz1)+ynorm(i,j-1,kz1))

             aixr = 0.5*(smrh(i+1,j,kz1)+smrh(i+1,j-1,kz1))*&
                    0.5*(mdot(i,j,kz1)+mdot(i,j-1,kz1))*&
                    0.5*(ynorm(i,j,kz1)+ynorm(i,j-1,kz1))

             aixl = 0.5*(smrh(i-1,j,kz1)+smrh(i-1,j-1,kz1))*&
                    0.5*(mdot(i,j,kz1)+mdot(i,j-1,kz1))*&
                    0.5*(ynorm(i,j,kz1)+ynorm(i,j-1,kz1))

             aiyr = 0.5*(smrh(i,j+1,kz1)+smrh(i,j,kz1))*&
                    0.5*(mdot(i,j,kz1)+mdot(i,j-1,kz1))*&
                    0.5*(ynorm(i,j,kz1)+ynorm(i,j-1,kz1))

             aiyl = 0.5*(smrh(i,j-1,kz1)+smrh(i,j-2,kz1))*&
                    0.5*(mdot(i,j,kz1)+mdot(i,j-1,kz1))*&
                    0.5*(ynorm(i,j,kz1)+ynorm(i,j-1,kz1))

             ! Method 1 - Smeared Density at Cell Faces

             !aicc = (rho1y(i,j,kz1)+rho1y(i,j,kz1))*&
             !       0.5*(mdot(i,j,kz1)+mdot(i,j-1,kz1))*&
             !       0.5*(ynorm(i,j,kz1)+ynorm(i,j-1,kz1))

             !aixr = (rho1y(i+1,j,kz1)+rho1y(i+1,j,kz1))*&
             !       0.5*(mdot(i,j,kz1)+mdot(i,j-1,kz1))*&
             !       0.5*(ynorm(i,j,kz1)+ynorm(i,j-1,kz1))

             !aixl = (rho1y(i-1,j,kz1)+rho1y(i-1,j,kz1))*&
             !       0.5*(mdot(i,j,kz1)+mdot(i,j-1,kz1))*&
             !       0.5*(ynorm(i,j,kz1)+ynorm(i,j-1,kz1))

             !aiyr = (rho1y(i,j+1,kz1)+rho1y(i,j+1,kz1))*&
             !       0.5*(mdot(i,j,kz1)+mdot(i,j-1,kz1))*&
             !       0.5*(ynorm(i,j,kz1)+ynorm(i,j-1,kz1))

             !aiyl = (rho1y(i,j-1,kz1)+rho1y(i,j-1,kz1))*&
             !       0.5*(mdot(i,j,kz1)+mdot(i,j-1,kz1))*&
             !       0.5*(ynorm(i,j,kz1)+ynorm(i,j-1,kz1))

             dsdxp = (aixr - aicc)*dx1
             dsdxm = (aicc - aixl)*dx1
             dsdyp = (aiyr - aicc)*dy1
             dsdym = (aicc - aiyl)*dy1

             !- Variable Viscosity Implementation (ru1 is 1/Re)
             txxp = ru1*0.25*(visc(i,j,kz1)+visc(i+1,j,kz1)+visc(i,j-1,kz1)+visc(i+1,j-1,kz1))*(dvdxp+dsdxp)
             txxm = ru1*0.25*(visc(i,j,kz1)+visc(i-1,j,kz1)+visc(i,j-1,kz1)+visc(i-1,j-1,kz1))*(dvdxm+dsdxm)
             tyyp = ru1*visc(i,j,kz1)  *(dvdyp+dsdyp)
             tyym = ru1*visc(i,j-1,kz1)*(dvdym+dsdym)

             Mdens = ( rho1y(i,j,kz1) + rho2y(i,j,kz1) )  ! Mixture inverse density.

             ! calculate RHS for v-momentum
             rv(i,j,kz1) = - (frx*uxplus - flx*uxminus)/dx &
                           - (fry*vyplus - fly*vyminus)/dy & 
                           + Mdens* (txxp - txxm)*dx1                 &! diffusion - normal terms
                           + Mdens* (tyyp - tyym)*dy1                 &
                           + gravY                                    ! kpd - gravity term                         
 
             !if(ycell .ge. Lb) rv(i,j,kz1) = rv(i,j,kz1) - ins_convvel(HIGH,JAXIS)*(vyplus - vyminus)*dy1

          enddo
       enddo

END SUBROUTINE ins_rhs2d_weno3

!-------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------!

SUBROUTINE ins_rhs3d_weno3(uni,vni,wni,tv,ru1,      &
                        ix1,ix2,jy1,jy2,kz1,kz2, &
                        dx,dy,dz,ru,rv,rw,visc,  &
                        rho1x,rho2x,rho1y,rho2y, &
                        rho1z,rho2z,gravX, gravY, gravZ,&
                        mdot,smrh,xnorm,ynorm,znorm,crv,temp,blockID)

  !*****************************************************************
  ! This subroutine computes the centered discretization of the RHS 
  ! of the momentum equation (advection + viscous terms) on a 
  ! staggered uniform grid based on the Paramesh grid structure.
  !
  ! Input:  uni,vni,wni = velocity at timestep n
  !         tv          = eddy viscosity
  !         ru1         = molecular viscosity !- kpd - Inverse Reynolds No
  !         ix1,ix2     = starting and ending x indices
  !         jy1,jy2     = starting and ending y indices
  !         kz1,kz2     = starting and ending z indices
  !         dx,dy,dz    = grid spacing in x, y, and z directions
  !
  ! Output: ru,rv,rw    = RHS of u, v, and w momentum equations
  !
  ! E. Balaras   July 1999
  ! P. Rabenold  August 2006
  !**************************************************************

      use Driver_interface, ONLY : Driver_abortFlash

      use RuntimeParameters_interface, ONLY : RuntimeParameters_get

      use IncompNS_data, ONLY : ins_iConvU,ins_convvel

      use Heat_AD_data, only: ht_Ra

      use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_getBlkCenterCoords, Grid_getDeltas

      use Simulation_data, only: sim_yMax, sim_yMin

      implicit none
      INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2, kz1, kz2
      REAL, INTENT(IN):: ru1, dx, dy, dz, gravX, gravY, gravZ
      REAL, DIMENSION(:,:,:), INTENT(IN):: uni, vni, wni, tv, visc, rho1x, rho2x, rho1y, rho2y 
      REAL, DIMENSION(:,:,:), INTENT(IN):: rho1z,rho2z 
      REAL, DIMENSION(:,:,:), INTENT(IN):: mdot,xnorm,ynorm,znorm,smrh,crv,temp
      REAL, DIMENSION(:,:,:), INTENT(OUT):: ru, rv, rw
      INTEGER, INTENT(IN) :: blockID

      INTEGER:: i, j, k
      REAL:: dx1, dy1, dz1, Mdens
      ! x-component variables
      REAL:: uxplus, uxminus, vxplus, vxminus, wxplus, wxminus, &
             uyplus, uyminus, uzplus, uzminus
      REAL:: dudxp, dudxm, dudyp, dudym, dudzp, dudzm, dvdxp, dvdxm, &
             dwdxp, dwdxm
      REAL:: tvjp, tvjm, tvkp, tvkm
      REAL:: txxp, txxm, tyyp, tyym, tzzp, tzzm
      REAL:: txyp, txym, txzp, txzm
      ! additional y-component variables
      REAL:: vyplus, vyminus, vzplus, vzminus, wyplus, wyminus
      REAL:: dvdyp, dvdym, dvdzp, dvdzm, dwdyp, dwdym
      REAL:: dsdxp, dsdxm, dsdyp, dsdym, dsdzp, dsdzm
      REAL:: aicc,aixr,aixl,aiyr,aiyl,aizr,aizl
      REAL:: tvip, tvim
      REAL:: tyzp, tyzm
      ! additional z-component variables
      REAL:: wzplus, wzminus
      REAL:: dwdzp, dwdzm, tempface

      REAL:: vvip,vvim,vvjp,vvjm,vvkp,vvkm

      real :: eps, &
              s1r,s2r,s3r,s4r,s5r,s1l,s2l,s3l,s4l,s5l, &
              rIS1r,rIS2r,rIS3r,rIS1l,rIS2l,rIS3l, &
              aT1r,aT2r,aT3r,aT1l,aT2l,aT3l, &
              a1r,a2r,a3r,a1l,a2l,a3l, &
              fT1r,fT2r,fT3r,fT1l,fT2l,fT3l, &
              frx,flx,fry,fly,frz,flz

      real :: ycell,Lb
      real :: del(MDIM),bsize(MDIM),coord(MDIM)
      real, dimension(2,MDIM) :: boundBox

      ! grid spacings
      dx1 = 1.0/dx
      dy1 = 1.0/dy
      dz1 = 1.0/dz

      eps = 1E-15

      Lb  = sim_yMax-sim_yMin-2.0

      call Grid_getDeltas(blockID,del)
      call Grid_getBlkCenterCoords(blockId,coord)
      call Grid_getBlkBoundBox(blockId,boundBox)

      bsize(:) = boundBox(2,:) - boundBox(1,:)

      !++++++++++  U-COMPONENT (Variable Density)  ++++++++++
      do k = kz1,kz2
         do j = jy1,jy2
            do i = ix1,ix2+1

              ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                       real(j - NGUARD - 1)*del(JAXIS)  +  &
                       0.5*del(JAXIS)
 
               ! Advection Terms

               uxplus  = (uni(i+1,j,k) + uni(i,j,k))*0.5 
               uxminus = (uni(i-1,j,k) + uni(i,j,k))*0.5

               vyplus  = (vni(i,j+1,k) + vni(i-1,j+1,k))*0.5
               vyminus = (vni(i,j,k) + vni(i-1,j,k))*0.5

               wzplus  = (wni(i,j,k+1) + wni(i-1,j,k+1))*0.5
               wzminus = (wni(i,j,k) + wni(i-1,j,k))*0.5

               uyplus = (uni(i,j+1,kz1) + uni(i,j,kz1))*0.5
               uyminus = (uni(i,j,kz1) + uni(i,j-1,kz1))*0.5

               !----------------- WENO3 X-Direction ------------!
               if (uxplus .gt. 0) then     ! u = (+) Downwind

                s1r = uni(i-2,j,k)
                s2r = uni(i-1,j,k)
                s3r = uni(i,j,k)
                s4r = uni(i+1,j,k)
                s5r = uni(i+2,j,k)

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

                s1r = uni(i-1,j,k)
                s2r = uni(i,j,k)
                s3r = uni(i+1,j,k)
                s4r = uni(i+2,j,k)
                s5r = uni(i+3,j,k)

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


               if (uxminus .gt. 0) then     ! u = (+) Downwind  

                s1l = uni(i-3,j,k)
                s2l = uni(i-2,j,k)
                s3l = uni(i-1,j,k)
                s4l = uni(i,j,k)
                s5l = uni(i+1,j,k)

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

                s1l = uni(i-2,j,k)
                s2l = uni(i-1,j,k)
                s3l = uni(i,j,k)
                s4l = uni(i+1,j,k)
                s5l = uni(i+2,j,k)

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
               !- WENO3 interpolated VEL values at cell center
               !---------------------------------------------------------
               frx = a1r*fT1r + a2r*fT2r + a3r*fT3r
               flx = a1l*fT1l + a2l*fT2l + a3l*fT3l
               !---------------------------------------------------------
               !---------------------------------------------------------

               !----------------- WENO3 Y-Direction ------------!
               if (vyplus .gt. 0) then     ! u = (+) Downwind

                s1r = uni(i,j-2,k)
                s2r = uni(i,j-1,k)
                s3r = uni(i,j,k)
                s4r = uni(i,j+1,k)
                s5r = uni(i,j+2,k)

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

                s1r = uni(i,j-1,k)
                s2r = uni(i,j,k)
                s3r = uni(i,j+1,k)
                s4r = uni(i,j+2,k)
                s5r = uni(i,j+3,k)

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

               if (vyminus .gt. 0) then     ! u = (+) Downwind

                s1l = uni(i,j-3,k)
                s2l = uni(i,j-2,k)
                s3l = uni(i,j-1,k)
                s4l = uni(i,j,k)
                s5l = uni(i,j+1,k)

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

                s1l = uni(i,j-2,k)
                s2l = uni(i,j-1,k)
                s3l = uni(i,j,k)
                s4l = uni(i,j+1,k)
                s5l = uni(i,j+2,k)

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
               !- WENO3 interpolated VEL values at cell center
               !---------------------------------------------------------
               fry = a1r*fT1r + a2r*fT2r + a3r*fT3r
               fly = a1l*fT1l + a2l*fT2l + a3l*fT3l
               !---------------------------------------------------------
               !---------------------------------------------------------

               !----------------- WENO3 Z-Direction ------------!
               if (wzplus .gt. 0) then     ! u = (+) Downwind

                s1r = uni(i,j,k-2)
                s2r = uni(i,j,k-1)
                s3r = uni(i,j,k)
                s4r = uni(i,j,k+1)
                s5r = uni(i,j,k+2)

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

                s1r = uni(i,j,k-1)
                s2r = uni(i,j,k)
                s3r = uni(i,j,k+1)
                s4r = uni(i,j,k+2)
                s5r = uni(i,j,k+3)

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

               if (wzminus .gt. 0) then     ! u = (+) Downwind

                s1l = uni(i,j,k-3)
                s2l = uni(i,j,k-2)
                s3l = uni(i,j,k-1)
                s4l = uni(i,j,k)
                s5l = uni(i,j,k+1)

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

                s1l = uni(i,j,k-2)
                s2l = uni(i,j,k-1)
                s3l = uni(i,j,k)
                s4l = uni(i,j,k+1)
                s5l = uni(i,j,k+2)

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
               !- WENO3 interpolated VEL values at cell center
               !---------------------------------------------------------
               frz = a1r*fT1r + a2r*fT2r + a3r*fT3r
               flz = a1l*fT1l + a2l*fT2l + a3l*fT3l
               !---------------------------------------------------------
               !---------------------------------------------------------

               ! Diffusion Terms

               ! get derivatives at 1/2 locations
               dudxp = (uni(i+1,j  ,k  ) - uni(i  ,j  ,k  ))*dx1
               dudxm = (uni(i  ,j  ,k  ) - uni(i-1,j  ,k  ))*dx1
               dudyp = (uni(i  ,j+1,k  ) - uni(i  ,j  ,k  ))*dy1
               dudym = (uni(i  ,j  ,k  ) - uni(i  ,j-1,k  ))*dy1
               dudzp = (uni(i  ,j  ,k+1) - uni(i  ,j  ,k  ))*dz1
               dudzm = (uni(i  ,j  ,k  ) - uni(i  ,j  ,k-1))*dz1
               dvdxp = (vni(i  ,j+1,k  ) - vni(i-1,j+1,k  ))*dx1
               dvdxm = (vni(i  ,j  ,k  ) - vni(i-1,j  ,k  ))*dx1
               dwdxp = (wni(i  ,j  ,k+1) - wni(i-1,j  ,k+1))*dx1
               dwdxm = (wni(i  ,j  ,k  ) - wni(i-1,j  ,k  ))*dx1

               ! velcoity jump source terms
               aicc = 0.5*(smrh(i,j,k)+smrh(i-1,j,k))*&
                      0.5*(mdot(i,j,k)+mdot(i-1,j,k))*&
                      0.5*(xnorm(i,j,k)+xnorm(i-1,j,k))

               aixr = 0.5*(smrh(i,j,k)+smrh(i+1,j,k))*&
                      0.5*(mdot(i,j,k)+mdot(i-1,j,k))*&
                      0.5*(xnorm(i,j,k)+xnorm(i-1,j,k))

               aixl = 0.5*(smrh(i-1,j,k)+smrh(i-2,j,k))*&
                      0.5*(mdot(i,j,k)+mdot(i-1,j,k))*&
                      0.5*(xnorm(i,j,k)+xnorm(i-1,j,k))

               aiyr = 0.5*(smrh(i,j+1,k)+smrh(i-1,j+1,k))*&
                      0.5*(mdot(i,j,k)+mdot(i-1,j,k))*&
                      0.5*(xnorm(i,j,k)+xnorm(i-1,j,k)) 

               aiyl = 0.5*(smrh(i,j-1,k)+smrh(i-1,j-1,k))*&
                      0.5*(mdot(i,j,k)+mdot(i-1,j,k))*&
                      0.5*(xnorm(i,j,k)+xnorm(i-1,j,k))
 
               aizr = 0.5*(smrh(i,j,k+1)+smrh(i-1,j,k+1))*&
                      0.5*(mdot(i,j,k)+mdot(i-1,j,k))*&
                      0.5*(xnorm(i,j,k)+xnorm(i-1,j,k)) 

               aizl = 0.5*(smrh(i,j,k-1)+smrh(i-1,j,k-1))*&
                      0.5*(mdot(i,j,k)+mdot(i-1,j,k))*&
                      0.5*(xnorm(i,j,k)+xnorm(i-1,j,k))
 
               dsdxp = (aixr - aicc)*dx1
               dsdxm = (aicc - aixl)*dx1
               dsdyp = (aiyr - aicc)*dy1
               dsdym = (aicc - aiyl)*dy1
               dsdzp = (aizr - aicc)*dz1
               dsdzm = (aicc - aizl)*dz1

               !- kpd - Eddy viscosity (Cell Center value) at diagonals
               tvjp = 0.25*(tv(i-1,j  ,k  ) + tv(i  ,j  ,k  ) + &
                            tv(i  ,j+1,k  ) + tv(i-1,j+1,k  ))
               tvjm = 0.25*(tv(i-1,j-1,k  ) + tv(i  ,j-1,k  ) + &
                            tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))
               tvkp = 0.25*(tv(i-1,j  ,k  ) + tv(i  ,j  ,k  ) + &
                            tv(i  ,j  ,k+1) + tv(i-1,j  ,k+1))
               tvkm = 0.25*(tv(i-1,j  ,k-1) + tv(i  ,j,  k-1) + &
                            tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))

               !- kpd - Molecular viscosity (Cell Center value) at diagonals
               vvip = visc(i,j,k)
               vvim = visc(i-1,j,k)
               vvjp = 0.25*(visc(i-1,j  ,k  ) + visc(i  ,j  ,k  ) + &
                            visc(i  ,j+1,k  ) + visc(i-1,j+1,k  ))
               vvjm = 0.25*(visc(i-1,j-1,k  ) + visc(i  ,j-1,k  ) + &
                            visc(i  ,j  ,k  ) + visc(i-1,j  ,k  ))
               vvkp = 0.25*(visc(i-1,j  ,k  ) + visc(i  ,j  ,k  ) + &
                            visc(i  ,j  ,k+1) + visc(i-1,j  ,k+1))
               vvkm = 0.25*(visc(i-1,j  ,k-1) + visc(i  ,j,  k-1) + &
                            visc(i  ,j  ,k  ) + visc(i-1,j  ,k  ))


               ! flux of normal total stresses, mu*dU/dXj (ru1 is 1/Re)
               txxp = (ru1*vvip + 2.0*tv(i,j,k))   *(dudxp+dsdxp)    ! KPD 2*nut, but not nu 
               txxm = (ru1*vvim + 2.0*tv(i-1,j,k)) *(dudxm+dsdxm)    ! KPD 2*nut, but not nu 
               tyyp = (ru1*vvjp + tvjp)            *(dudyp+dsdyp)
               tyym = (ru1*vvjm + tvjm)            *(dudym+dsdym)
               tzzp = (ru1*vvkp + tvkp)            *(dudzp+dsdzp)
               tzzm = (ru1*vvkm + tvkm)            *(dudzm+dsdzm)

               ! flux of cross SGS stresses, mu*dUi/dX
               txyp = tvjp*dvdxp                             ! KPD no mol visc, turb only 
               txym = tvjm*dvdxm                             ! KPD no mol visc, turb only 
               txzp = tvkp*dwdxp                             ! KPD no mol visc, turb only 
               txzm = tvkm*dwdxm                             ! KPD no mol visc, turb only 

               !- kpd - Mixture inverse density
               Mdens = ( rho1x(i,j,k) + rho2x(i,j,k) )

               !- AD - temperature at face for natural convection
               tempface = 0.5*(temp(i,j,k) + temp(i-1,j,k))

               ! calculate RHS for u-momentum
               ru(i,j,k) = - (frx*uxplus - flx*uxminus)/dx &
                           - (fry*vyplus - fly*vyminus)/dy &
                           - (frz*wzplus - flz*wzminus)/dz &
                           + Mdens*(txxp - txxm)*dx1                &! diffusion - normal terms
                           + Mdens*(tyyp - tyym)*dy1                &
                           + Mdens*(tzzp - tzzm)*dz1                &
                           + Mdens*(txyp - txym)*dy1                &! TURBULENT cross terms
                           + Mdens*(txzp - txzm)*dz1                &
                           + gravX*(1.0  + ht_Ra*tempface)  
 
               
               if(ycell .ge. Lb) ru(i,j,k) = ru(i,j,k) - ins_convvel(HIGH,JAXIS)*(uyplus - uyminus)*dy1

            enddo
         enddo
      enddo

      !++++++++++  V-COMPONENT  ++++++++++

      do k = kz1,kz2
         do j = jy1,jy2+1
            do i = ix1,ix2

               ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                       real(j - NGUARD - 1)*del(JAXIS)

              ! Advection Terms

               uxplus  = (uni(i+1,j,k) + uni(i+1,j-1,k))*0.5 
               uxminus = (uni(i,j,k) + uni(i,j-1,k))*0.5

               vyplus  = (vni(i,j+1,k) + vni(i,j,k))*0.5
               vyminus = (vni(i,j-1,k) + vni(i,j,k))*0.5

               wzplus  = (wni(i,j,k+1) + wni(i,j-1,k+1))*0.5
               wzminus = (wni(i,j,k) + wni(i,j-1,k))*0.5

               !----------------- WENO3 X-Direction ------------!
               if (uxplus .gt. 0) then     ! u = (+) Downwind

                s1r = vni(i-2,j,k)
                s2r = vni(i-1,j,k)
                s3r = vni(i,j,k)
                s4r = vni(i+1,j,k)
                s5r = vni(i+2,j,k)

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

                s1r = vni(i-1,j,k)
                s2r = vni(i,j,k)
                s3r = vni(i+1,j,k)
                s4r = vni(i+2,j,k)
                s5r = vni(i+3,j,k)

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


               if (uxminus .gt. 0) then     ! u = (+) Downwind  

                s1l = vni(i-3,j,k)
                s2l = vni(i-2,j,k)
                s3l = vni(i-1,j,k)
                s4l = vni(i,j,k)
                s5l = vni(i+1,j,k)

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

                s1l = vni(i-2,j,k)
                s2l = vni(i-1,j,k)
                s3l = vni(i,j,k)
                s4l = vni(i+1,j,k)
                s5l = vni(i+2,j,k)

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
               !- WENO3 interpolated VEL values at cell center
               !---------------------------------------------------------
               frx = a1r*fT1r + a2r*fT2r + a3r*fT3r
               flx = a1l*fT1l + a2l*fT2l + a3l*fT3l
               !---------------------------------------------------------
               !---------------------------------------------------------

               !----------------- WENO3 Y-Direction ------------!
               if (vyplus .gt. 0) then     ! u = (+) Downwind

                s1r = vni(i,j-2,k)
                s2r = vni(i,j-1,k)
                s3r = vni(i,j,k)
                s4r = vni(i,j+1,k)
                s5r = vni(i,j+2,k)

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

                s1r = vni(i,j-1,k)
                s2r = vni(i,j,k)
                s3r = vni(i,j+1,k)
                s4r = vni(i,j+2,k)
                s5r = vni(i,j+3,k)

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

               if (vyminus .gt. 0) then     ! u = (+) Downwind

                s1l = vni(i,j-3,k)
                s2l = vni(i,j-2,k)
                s3l = vni(i,j-1,k)
                s4l = vni(i,j,k)
                s5l = vni(i,j+1,k)

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

                s1l = vni(i,j-2,k)
                s2l = vni(i,j-1,k)
                s3l = vni(i,j,k)
                s4l = vni(i,j+1,k)
                s5l = vni(i,j+2,k)

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
               !- WENO3 interpolated VEL values at cell center
               !---------------------------------------------------------
               fry = a1r*fT1r + a2r*fT2r + a3r*fT3r
               fly = a1l*fT1l + a2l*fT2l + a3l*fT3l
               !---------------------------------------------------------
               !---------------------------------------------------------

               !----------------- WENO3 Z-Direction ------------!
               if (wzplus .gt. 0) then     ! u = (+) Downwind

                s1r = vni(i,j,k-2)
                s2r = vni(i,j,k-1)
                s3r = vni(i,j,k)
                s4r = vni(i,j,k+1)
                s5r = vni(i,j,k+2)

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

                s1r = vni(i,j,k-1)
                s2r = vni(i,j,k)
                s3r = vni(i,j,k+1)
                s4r = vni(i,j,k+2)
                s5r = vni(i,j,k+3)

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

               if (wzminus .gt. 0) then     ! u = (+) Downwind

                s1l = vni(i,j,k-3)
                s2l = vni(i,j,k-2)
                s3l = vni(i,j,k-1)
                s4l = vni(i,j,k)
                s5l = vni(i,j,k+1)

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

                s1l = vni(i,j,k-2)
                s2l = vni(i,j,k-1)
                s3l = vni(i,j,k)
                s4l = vni(i,j,k+1)
                s5l = vni(i,j,k+2)

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
               !- WENO3 interpolated VEL values at cell center
               !---------------------------------------------------------
               frz = a1r*fT1r + a2r*fT2r + a3r*fT3r
               flz = a1l*fT1l + a2l*fT2l + a3l*fT3l
               !---------------------------------------------------------
               !---------------------------------------------------------

               ! Diffusion Terms

               ! get derivatives at 1/2 locations
               dvdxp = (vni(i+1,j  ,k  ) - vni(i  ,j  ,k  ))*dx1
               dvdxm = (vni(i  ,j  ,k  ) - vni(i-1,j  ,k  ))*dx1
               dvdyp = (vni(i  ,j+1,k  ) - vni(i  ,j  ,k  ))*dy1
               dvdym = (vni(i  ,j  ,k  ) - vni(i  ,j-1,k  ))*dy1
               dvdzp = (vni(i  ,j  ,k+1) - vni(i  ,j  ,k  ))*dz1
               dvdzm = (vni(i  ,j  ,k  ) - vni(i  ,j  ,k-1))*dz1
               dudyp = (uni(i+1,j  ,k  ) - uni(i+1,j-1,k  ))*dy1
               dudym = (uni(i  ,j  ,k  ) - uni(i  ,j-1,k  ))*dy1
               dwdyp = (wni(i  ,j  ,k+1) - wni(i  ,j-1,k+1))*dy1
               dwdym = (wni(i  ,j  ,k  ) - wni(i  ,j-1,k  ))*dy1

               ! velcoity jump source terms

               aicc = 0.5*(smrh(i,j,k)+smrh(i,j-1,k))*&
                      0.5*(mdot(i,j,k)+mdot(i,j-1,k))*&
                      0.5*(ynorm(i,j,k)+ynorm(i,j-1,k))

               aixr = 0.5*(smrh(i+1,j,k)+smrh(i+1,j-1,k))*&
                      0.5*(mdot(i,j,k)+mdot(i,j-1,k))*&
                      0.5*(ynorm(i,j,k)+ynorm(i,j-1,k))

               aixl = 0.5*(smrh(i-1,j,k)+smrh(i-1,j-1,k))*&
                      0.5*(mdot(i,j,k)+mdot(i,j-1,k))*&
                      0.5*(ynorm(i,j,k)+ynorm(i,j-1,k))

               aiyr = 0.5*(smrh(i,j+1,k)+smrh(i,j,k))*&
                      0.5*(mdot(i,j,k)+mdot(i,j-1,k))*&
                      0.5*(ynorm(i,j,k)+ynorm(i,j-1,k))

               aiyl = 0.5*(smrh(i,j-1,k)+smrh(i,j-2,k))*&
                      0.5*(mdot(i,j,k)+mdot(i,j-1,k))*&
                      0.5*(ynorm(i,j,k)+ynorm(i,j-1,k))

               aizr = 0.5*(smrh(i,j,k+1)+smrh(i,j-1,k+1))*&
                      0.5*(mdot(i,j,k)+mdot(i,j-1,k))*&
                      0.5*(ynorm(i,j,k)+ynorm(i,j-1,k))

               aizl = 0.5*(smrh(i,j,k-1)+smrh(i,j-1,k-1))*&
                      0.5*(mdot(i,j,k)+mdot(i,j-1,k))*&
                      0.5*(ynorm(i,j,k)+ynorm(i,j-1,k))

               dsdxp = (aixr - aicc)*dx1
               dsdxm = (aicc - aixl)*dx1
               dsdyp = (aiyr - aicc)*dy1
               dsdym = (aicc - aiyl)*dy1
               dsdzp = (aizr - aicc)*dz1
               dsdzm = (aicc - aizl)*dz1

               !- kpd - Eddy viscosity (Cell Center value) at diagonals
               tvip = 0.25*(tv(i  ,j-1,k  ) + tv(i+1,j-1,k  ) + &
                            tv(i+1,j  ,k  ) + tv(i  ,j  ,k  ))
               tvim = 0.25*(tv(i-1,j-1,k  ) + tv(i  ,j-1,k  ) + &
                            tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))
               tvkp = 0.25*(tv(i  ,j-1,k  ) + tv(i  ,j-1,k+1) + &
                            tv(i  ,j  ,k+1) + tv(i  ,j  ,k  ))
               tvkm = 0.25*(tv(i  ,j-1,k-1) + tv(i  ,j-1,k  ) + &
                            tv(i  ,j  ,k  ) + tv(i  ,j  ,k-1))

               !- kpd - Molecular viscosity (Cell Center value) at diagonals
               vvip = 0.25*(visc(i  ,j-1,k  ) + visc(i+1,j-1,k  ) + &
                            visc(i+1,j  ,k  ) + visc(i  ,j  ,k  ))
               vvim = 0.25*(visc(i-1,j-1,k  ) + visc(i  ,j-1,k  ) + &
                            visc(i  ,j  ,k  ) + visc(i-1,j  ,k  ))
               vvjp = visc(i,j,k) 
               vvjm = visc(i,j-1,k)
               vvkp = 0.25*(visc(i  ,j-1,k  ) + visc(i  ,j-1,k+1) + &
                            visc(i  ,j  ,k+1) + visc(i  ,j  ,k  ))
               vvkm = 0.25*(visc(i  ,j-1,k-1) + visc(i  ,j-1,k  ) + &
                            visc(i  ,j  ,k  ) + visc(i  ,j  ,k-1))

               ! flux of normal total stresses (ru1 is 1/Re)
               txxp = (ru1*vvip + tvip           )*(dvdxp+dsdxp)
               txxm = (ru1*vvim + tvim           )*(dvdxm+dsdxm)
               tyyp = (ru1*vvjp + 2.0*tv(i,j,k)  )*(dvdyp+dsdyp)
               tyym = (ru1*vvjm + 2.0*tv(i,j-1,k))*(dvdym+dsdym)
               tzzp = (ru1*vvkp + tvkp           )*(dvdzp+dsdzp)
               tzzm = (ru1*vvkm + tvkm           )*(dvdzm+dsdzm)

               ! flux of cross SGS stresses
               txyp = tvip*dudyp
               txym = tvim*dudym
               tyzp = tvkp*dwdyp
               tyzm = tvkm*dwdym

               !- kpd - Mixture inverse density
               Mdens = ( rho1y(i,j,k) + rho2y(i,j,k) ) 

               !- AD - temperature at face for natural convection
               tempface = 0.5*(temp(i,j,k) + temp(i,j-1,k))

               ! calculate RHS for v-momentum
               rv(i,j,k) = - (frx*uxplus - flx*uxminus)/dx &
                           - (fry*vyplus - fly*vyminus)/dy &
                           - (frz*wzplus - flz*wzminus)/dz &
                           + Mdens* (txxp - txxm)*dx1                &! diffusion - normal terms
                           + Mdens* (tyyp - tyym)*dy1                &
                           + Mdens* (tzzp - tzzm)*dz1                &
                           + Mdens* (txyp - txym)*dx1                &! diffusion - cross terms
                           + Mdens* (tyzp - tyzm)*dz1                &
                           + gravY*(1.0  + ht_Ra*tempface)           

              if(ycell .ge. Lb) rv(i,j,k) = rv(i,j,k) - ins_convvel(HIGH,JAXIS)*(vyplus - vyminus)*dy1

            enddo
         enddo
      enddo

      !++++++++++  W-COMPONENT  ++++++++++
      
      do k = kz1,kz2+1
         do j = jy1,jy2
            do i = ix1,ix2

              ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                       real(j - NGUARD - 1)*del(JAXIS)  +  &
                       0.5*del(JAXIS)
 
               ! Advection Terms

               uxplus  = (uni(i+1,j,k) + uni(i+1,j,k-1))*0.5 
               uxminus = (uni(i,j,k) + uni(i,j,k-1))*0.5

               vyplus  = (vni(i,j+1,k) + vni(i,j+1,k-1))*0.5
               vyminus = (vni(i,j,k) + vni(i,j,k-1))*0.5

               wzplus  = (wni(i,j,k+1) + wni(i,j,k))*0.5
               wzminus = (wni(i,j,k-1) + wni(i,j,k))*0.5

               wyplus  = (wni(i,j+1,kz1) + wni(i,j,kz1))*0.5
               wyminus = (wni(i,j,kz1)   + wni(i,j-1,kz1))*0.5

               !----------------- WENO3 X-Direction ------------!
               if (uxplus .gt. 0) then     ! u = (+) Downwind

                s1r = wni(i-2,j,k)
                s2r = wni(i-1,j,k)
                s3r = wni(i,j,k)
                s4r = wni(i+1,j,k)
                s5r = wni(i+2,j,k)

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

                s1r = wni(i-1,j,k)
                s2r = wni(i,j,k)
                s3r = wni(i+1,j,k)
                s4r = wni(i+2,j,k)
                s5r = wni(i+3,j,k)

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


               if (uxminus .gt. 0) then     ! u = (+) Downwind  

                s1l = wni(i-3,j,k)
                s2l = wni(i-2,j,k)
                s3l = wni(i-1,j,k)
                s4l = wni(i,j,k)
                s5l = wni(i+1,j,k)

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

                s1l = wni(i-2,j,k)
                s2l = wni(i-1,j,k)
                s3l = wni(i,j,k)
                s4l = wni(i+1,j,k)
                s5l = wni(i+2,j,k)

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
               !- WENO3 interpolated VEL values at cell center
               !---------------------------------------------------------
               frx = a1r*fT1r + a2r*fT2r + a3r*fT3r
               flx = a1l*fT1l + a2l*fT2l + a3l*fT3l
               !---------------------------------------------------------
               !---------------------------------------------------------

               !----------------- WENO3 Y-Direction ------------!
               if (vyplus .gt. 0) then     ! u = (+) Downwind

                s1r = wni(i,j-2,k)
                s2r = wni(i,j-1,k)
                s3r = wni(i,j,k)
                s4r = wni(i,j+1,k)
                s5r = wni(i,j+2,k)

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

                s1r = wni(i,j-1,k)
                s2r = wni(i,j,k)
                s3r = wni(i,j+1,k)
                s4r = wni(i,j+2,k)
                s5r = wni(i,j+3,k)

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

               if (vyminus .gt. 0) then     ! u = (+) Downwind

                s1l = wni(i,j-3,k)
                s2l = wni(i,j-2,k)
                s3l = wni(i,j-1,k)
                s4l = wni(i,j,k)
                s5l = wni(i,j+1,k)

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

                s1l = wni(i,j-2,k)
                s2l = wni(i,j-1,k)
                s3l = wni(i,j,k)
                s4l = wni(i,j+1,k)
                s5l = wni(i,j+2,k)

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
               !- WENO3 interpolated VEL values at cell center
               !---------------------------------------------------------
               fry = a1r*fT1r + a2r*fT2r + a3r*fT3r
               fly = a1l*fT1l + a2l*fT2l + a3l*fT3l
               !---------------------------------------------------------
               !---------------------------------------------------------

               !----------------- WENO3 Z-Direction ------------!
               if (wzplus .gt. 0) then     ! u = (+) Downwind

                s1r = wni(i,j,k-2)
                s2r = wni(i,j,k-1)
                s3r = wni(i,j,k)
                s4r = wni(i,j,k+1)
                s5r = wni(i,j,k+2)

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

                s1r = wni(i,j,k-1)
                s2r = wni(i,j,k)
                s3r = wni(i,j,k+1)
                s4r = wni(i,j,k+2)
                s5r = wni(i,j,k+3)

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

               if (wzminus .gt. 0) then     ! u = (+) Downwind

                s1l = wni(i,j,k-3)
                s2l = wni(i,j,k-2)
                s3l = wni(i,j,k-1)
                s4l = wni(i,j,k)
                s5l = wni(i,j,k+1)

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

                s1l = wni(i,j,k-2)
                s2l = wni(i,j,k-1)
                s3l = wni(i,j,k)
                s4l = wni(i,j,k+1)
                s5l = wni(i,j,k+2)

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
               !- WENO3 interpolated VEL values at cell center
               !---------------------------------------------------------
               frz = a1r*fT1r + a2r*fT2r + a3r*fT3r
               flz = a1l*fT1l + a2l*fT2l + a3l*fT3l
               !---------------------------------------------------------
               !---------------------------------------------------------

               ! Diffusion Terms

               ! get derivatives at 1/2 locations
               dwdxp = (wni(i+1,j  ,k  ) - wni(i  ,j  ,k  ))*dx1
               dwdxm = (wni(i  ,j  ,k  ) - wni(i-1,j  ,k  ))*dx1
               dwdyp = (wni(i  ,j+1,k  ) - wni(i  ,j  ,k  ))*dy1
               dwdym = (wni(i  ,j  ,k  ) - wni(i  ,j-1,k  ))*dy1
               dwdzp = (wni(i  ,j  ,k+1) - wni(i  ,j  ,k  ))*dz1
               dwdzm = (wni(i  ,j  ,k  ) - wni(i  ,j  ,k-1))*dz1
               dudzp = (uni(i+1,j  ,k  ) - uni(i+1,j  ,k-1))*dz1
               dudzm = (uni(i  ,j  ,k  ) - uni(i  ,j  ,k-1))*dz1
               dvdzp = (vni(i  ,j+1,k  ) - vni(i  ,j+1,k-1))*dz1
               dvdzm = (vni(i  ,j  ,k  ) - vni(i  ,j  ,k-1))*dz1

               ! velcoity jump source terms
               aicc = 0.5*(smrh(i,j,k)+smrh(i,j,k-1))*&
                      0.5*(mdot(i,j,k)+mdot(i,j,k-1))*&
                      0.5*(znorm(i,j,k)+znorm(i,j,k-1))

               aixr = 0.5*(smrh(i+1,j,k)+smrh(i+1,j,k-1))*&
                      0.5*(mdot(i,j,k)+mdot(i,j,k-1))*&
                      0.5*(znorm(i,j,k)+znorm(i,j,k-1))

               aixl = 0.5*(smrh(i-1,j,k)+smrh(i-1,j,k-1))*&
                      0.5*(mdot(i,j,k)+mdot(i,j,k-1))*&
                      0.5*(znorm(i,j,k)+znorm(i,j,k-1))

               aiyr = 0.5*(smrh(i,j+1,k)+smrh(i,j+1,k-1))*&
                      0.5*(mdot(i,j,k)+mdot(i,j,k-1))*&
                      0.5*(znorm(i,j,k)+znorm(i,j,k-1)) 

               aiyl = 0.5*(smrh(i,j-1,k)+smrh(i,j-1,k-1))*&
                      0.5*(mdot(i,j,k)+mdot(i,j,k-1))*&
                      0.5*(znorm(i,j,k)+znorm(i,j,k-1))
 
               aizr = 0.5*(smrh(i,j,k+1)+smrh(i,j,k))*&
                      0.5*(mdot(i,j,k)+mdot(i,j,k-1))*&
                      0.5*(znorm(i,j,k)+znorm(i,j,k-1))

               aizl = 0.5*(smrh(i,j,k-1)+smrh(i,j,k-2))*&
                      0.5*(mdot(i,j,k)+mdot(i,j,k-1))*&
                      0.5*(znorm(i,j,k)+znorm(i,j,k-1))

               dsdxp = (aixr - aicc)*dx1
               dsdxm = (aicc - aixl)*dx1
               dsdyp = (aiyr - aicc)*dy1
               dsdym = (aicc - aiyl)*dy1
               dsdzp = (aizr - aicc)*dz1
               dsdzm = (aicc - aizl)*dz1

               !- kpd - Eddy viscosity (Cell Center value) at diagonals
               tvip = 0.25*(tv(i  ,j  ,k-1) + tv(i+1,j  ,k-1) + &
                            tv(i+1,j  ,k  ) + tv(i  ,j  ,k  ))
               tvim = 0.25*(tv(i-1,j  ,k-1) + tv(i  ,j  ,k-1) + &
                            tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))
               tvjp = 0.25*(tv(i  ,j  ,k-1) + tv(i  ,j  ,k  ) + &
                            tv(i  ,j+1,k  ) + tv(i  ,j+1,k-1))
               tvjm = 0.25*(tv(i  ,j-1,k-1) + tv(i  ,j-1,k  ) + & 
                            tv(i  ,j  ,k  ) + tv(i  ,j  ,k-1))

               !- kpd - Molecular viscosity (Cell Center value) at diagonals
               vvip = 0.25*(visc(i  ,j  ,k-1) + visc(i+1,j  ,k-1) + &
                            visc(i+1,j  ,k  ) + visc(i  ,j  ,k  ))
               vvim = 0.25*(visc(i-1,j  ,k-1) + visc(i  ,j  ,k-1) + &
                            visc(i  ,j  ,k  ) + visc(i-1,j  ,k  ))
               vvjp = 0.25*(visc(i  ,j  ,k-1) + visc(i  ,j  ,k  ) + &
                            visc(i  ,j+1,k  ) + visc(i  ,j+1,k-1))
               vvjm = 0.25*(visc(i  ,j-1,k-1) + visc(i  ,j-1,k  ) + &
                            visc(i  ,j  ,k  ) + visc(i  ,j  ,k-1))
               vvkp = visc(i,j,k)
               vvkm = visc(i,j,k-1)

               ! flux of normal total stresses (ru1 is 1/Re)
               txxp = (ru1*vvip + tvip           )*(dwdxp+dsdxp)
               txxm = (ru1*vvim + tvim           )*(dwdxm+dsdxm)
               tyyp = (ru1*vvjp + tvjp           )*(dwdyp+dsdyp)
               tyym = (ru1*vvjm + tvjm           )*(dwdym+dsdym)
               tzzp = (ru1*vvkp + 2.0*tv(i,j,k)  )*(dwdzp+dsdzp)
               tzzm = (ru1*vvkm + 2.0*tv(i,j,k-1))*(dwdzm+dsdzm)

               ! flux of cross SGS stresses
               txzp = tvip*dudzp
               txzm = tvim*dudzm
               tyzp = tvjp*dvdzp
               tyzm = tvjm*dvdzm

               !- kpd - Mixture inverse density
               Mdens = ( rho1z(i,j,k) + rho2z(i,j,k) )

               !- AD - temperature at face for natural convection
               tempface = 0.5*(temp(i,j,k) + temp(i,j,k-1))

               ! calculate RHS for w-momentum
               rw(i,j,k) = - (frx*uxplus - flx*uxminus)/dx &
                           - (fry*vyplus - fly*vyminus)/dy &
                           - (frz*wzplus - flz*wzminus)/dz &
                           + Mdens* (txxp - txxm)*dx1                &! diffusion - normal terms
                           + Mdens* (tyyp - tyym)*dy1                &
                           + Mdens* (tzzp - tzzm)*dz1                &
                           + Mdens* (txzp - txzm)*dx1                &! diffusion - cross terms
                           + Mdens* (tyzp - tyzm)*dy1                &
                           + gravZ* (1.0  + ht_Ra*tempface)                  

              if(ycell .ge. Lb) rw(i,j,k) = rw(i,j,k) - ins_convvel(HIGH,JAXIS)*(wyplus - wyminus)*dy1
 
            enddo
         enddo
      enddo

END SUBROUTINE ins_rhs3d_weno3
