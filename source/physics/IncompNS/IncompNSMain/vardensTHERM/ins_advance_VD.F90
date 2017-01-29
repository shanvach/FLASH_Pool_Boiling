SUBROUTINE ins_predictor_VD(uni,vni,wni,unew,vnew,wnew,uold,vold,&
        wold,p,dt,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,gama,rhoa,alfa)

      ! This routine computes the intermediate velocities based on
      ! the explicit second-order Adams-Bashforth scheme (gamma=1.5,
      ! rhoa=-0.5,alfa=1), or a third order Runge-Kutta method.

      use IncompNS_data, ONLY : ins_prescoeff

      implicit none

#include "Flash.h"

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dt,dx,dy,dz
      REAL, DIMENSION(:,:,:), INTENT(IN) :: unew,vnew,wnew,&
                                            uold,vold,wold,&
                                            p
      REAL, DIMENSION(:,:,:), INTENT(IN OUT) :: uni,vni,wni

      REAL :: gama,rhoa,alfa

      !- kpd - Output of pressure solution method to the screen
      !if (ins_prescoeff .gt. 0.) then
      !   print*,"KPD - Pressure Correction scheme IS in use (ins_advance.F90)."
      !else
      !   print*,"KPD - Pressure Correction scheme is NOT in use (ins_advance.F90)."
      !end if

      !---------------------------------------------------------------------
      !---------------------------------------------------------------------
      !- kpd - When pressure correction scheme is not used ins_prescoeff = 0
      !           and dP/dx is not used in u* predictor step
      !        uni and vni are u* and v*
      !---------------------------------------------------------------------
      !---------------------------------------------------------------------

      uni(ix1:ix2+1,jy1:jy2,kz1:kz2) = uni(ix1:ix2+1,jy1:jy2,kz1:kz2) + &
         dt*(gama*unew(ix1:ix2+1,jy1:jy2,kz1:kz2) +                     &
             rhoa*uold(ix1:ix2+1,jy1:jy2,kz1:kz2) -                     &
             ins_prescoeff*alfa*( p(ix1:ix2+1,jy1:jy2,kz1:kz2) -        &
                                  p(ix1-1:ix2,jy1:jy2,kz1:kz2) )/dx)

      vni(ix1:ix2,jy1:jy2+1,kz1:kz2) = vni(ix1:ix2,jy1:jy2+1,kz1:kz2) + &
         dt*(gama*vnew(ix1:ix2,jy1:jy2+1,kz1:kz2) +                     &
             rhoa*vold(ix1:ix2,jy1:jy2+1,kz1:kz2) -                     &
             ins_prescoeff*alfa*( p(ix1:ix2,jy1:jy2+1,kz1:kz2) -        &
                                  p(ix1:ix2,jy1-1:jy2,kz1:kz2) )/dy)

#if NDIM == 3
      wni(ix1:ix2,jy1:jy2,kz1:kz2+1) = wni(ix1:ix2,jy1:jy2,kz1:kz2+1) + &
         dt*(gama*wnew(ix1:ix2,jy1:jy2,kz1:kz2+1) +                     &
             rhoa*wold(ix1:ix2,jy1:jy2,kz1:kz2+1) -                     &
             ins_prescoeff*alfa*( p(ix1:ix2,jy1:jy2,kz1:kz2+1) -        &
                                  p(ix1:ix2,jy1:jy2,kz1-1:kz2) )/dz)
#endif


END SUBROUTINE ins_predictor_VD

!########################################################################

SUBROUTINE ins_divergence_VD(uni,vni,wni,ix1,ix2,jy1,jy2,kz1,kz2,&
         dx,dy,dz,divv)

      ! This routine computes the divergence of the velocity field.

      implicit none

#include "Flash.h"

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dx,dy,dz
      REAL, DIMENSION(:,:,:), INTENT(IN) :: uni,vni,wni
      REAL, DIMENSION(:,:,:), INTENT(OUT) :: divv


      divv(ix1:ix2,jy1:jy2,kz1:kz2) =           & 
         ( uni(ix1+1:ix2+1,jy1:jy2,kz1:kz2) -   &
           uni(ix1:ix2,jy1:jy2,kz1:kz2) )/dx +  &       
         ( vni(ix1:ix2,jy1+1:jy2+1,kz1:kz2) -   &
           vni(ix1:ix2,jy1:jy2,kz1:kz2) )/dy


#if NDIM == 3
      divv(ix1:ix2,jy1:jy2,kz1:kz2) =           & 
           divv(ix1:ix2,jy1:jy2,kz1:kz2)     +  &
         ( wni(ix1:ix2,jy1:jy2,kz1+1:kz2+1) -   &
           wni(ix1:ix2,jy1:jy2,kz1:kz2) )/dz 
#endif

END SUBROUTINE ins_divergence_VD

!########################################################################

SUBROUTINE ins_divergence_VD_ML(uni,vni,wni,ix1,ix2,jy1,jy2,kz1,kz2,&
                                ix1gc,ix2gc,jy1gc,jy2gc,kz1gc,kz2gc,&
                                       dx,dy,dz,divv,s,pf,xnorm,ynorm,smrh,mdot,rho1,rho2,testval,rho1x,rho2x,rho1y,rho2y)

      ! This routine computes the divergence of the velocity field.

      implicit none

#include "Flash.h"

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dx,dy,dz
      REAL, DIMENSION(:,:,:), INTENT(IN) :: uni,vni,wni
      REAL, DIMENSION(:,:,:), INTENT(OUT) :: divv

! added by mslee
      integer :: i,j,k
      real :: sp,rhoxr,rhoxl,rhoyr,rhoyl,aixr,aixl,aiyr,aiyl
      INTEGER, INTENT(IN) :: ix1gc,ix2gc,jy1gc,jy2gc,kz1gc,kz2gc
      real, dimension(:,:,:), intent(in) :: s,pf
      real, intent(in) :: rho1,rho2
      real, dimension(:,:,:), intent(in) :: xnorm,ynorm,smrh,mdot,rho1x,rho2x,rho1y,rho2y
      real, dimension(:,:,:), intent(out) :: testval
!      real, dimension(ix1gc:ix2gc,jy1gc:jy2gc,kz1gc:kz2gc) :: smrh
!      real, allocatable :: xnorm(:,:,:),ynorm(:,:,:),rhof(:,:,:)
!      allocate(xnorm(ix1gc:ix2gc,jy1gc:jy2gc,kz1gc:kz2gc),ynorm(ix1gc:ix2gc,jy1gc:jy2gc,kz1gc:kz2gc),rhof(ix1gc:ix2gc,jy1gc:jy2gc,kz1gc:kz2gc))
!      real, dimension(:,:,:) :: pf,testsourcex,testsourcey


!         k=1
!         do j = jy1-1,jy2+1
!            do i = ix1-1,ix2+1
!               pf(i,j,k) = 0.
! 
!               if(s(i,j,k).ge.0.) then
!                  pf(i,j,k) = 1.                       
!               end if
!            end do
!         end do
!
!         k=1
!         do j = jy1-1,jy2+1
!            do i = ix1-1,ix2+1
!                if(pf(i,j,k).eq.0..and.pf(i,j+1,k).eq.1.) then
!                   testsourcey(i,j,k) = 0.1
!                end if
!                if(pf(i,j,k).eq.0..and.pf(i+1,j,k).eq.1.) then
!                   testsourcex(i,j,k) = 0.1
!                end if
!            end do
!         end do

!        sxl = s(i-1,j,k)
!        sxr = s(i+1,j,k)
!        syl = s(i,j-1,k)
!        syr = s(i,j+1,k)
!
!        adf = sqrt( ((sxr-sxl)/2./dx)**2 + ((syr-syl)/2./dy)**2 )
!        xnorm(i,j,k) = (sxr-sxl)/2./dx/adf
!        ynorm(i,j,k) = (syr-syl)/2./dy/adf

!!        xnorm(ix1:ix2,jy1:jy2,kz1:kz2) =           &
!!          (( s(ix1+1:ix2+1,jy1:jy2,kz1:kz2) -   &
!!            s(ix1-1:ix2-1,jy1:jy2,kz1:kz2) )/2./dx)/ &
!!            sqrt( ((s(ix1+1:ix2+1,jy1:jy2,kz1:kz2) - &
!!            s(ix1-1:ix2-1,jy1:jy2,kz1:kz2))/2./dx)**2 &
!!            + ((s(ix1:ix2,jy1+1:jy2+1,kz1:kz2) - &
!!            s(ix1:ix2,jy1-1:jy2-1,kz1:kz2))/2./dy)**2 )
!!            
!!        ynorm(ix1:ix2,jy1:jy2,kz1:kz2) =           &
!!          (( s(ix1:ix2,jy1+1:jy2+1,kz1:kz2) -   &
!!            s(ix1:ix2,jy1-1:jy2-1,kz1:kz2) )/2./dy)/ &
!!            sqrt( ((s(ix1+1:ix2+1,jy1:jy2,kz1:kz2) - &
!!            s(ix1-1:ix2-1,jy1:jy2,kz1:kz2))/2./dx)**2 &
!!            + ((s(ix1:ix2,jy1+1:jy2+1,kz1:kz2) - &
!!            s(ix1:ix2,jy1-1:jy2-1,kz1:kz2))/2./dy)**2 )
!!
!!        sp=1.d0*min(dx,dy)
!!        smhv(ix1:ix2,jy1:jy2,kz1:kz2) = &
!!                (1.d0 + derf(s(ix1:ix2,jy1:jy2,kz1:kz2)/sp))/2.d0

!        smrh(ix1:ix2,jy1:jy2,kz1:kz2) = &
!                rho2 + (rho1 - rho2) * smhv(ix1:ix2,jy1:jy2,kz1:kz2)

!                1./rho1 + pf(ix1:ix2,jy1:jy2,kz1:kz2)*(1./rho2 - 1./rho1)


!           smhv(ix1:ix2,jy1:jy2,kz1:kz2) = &
!                1./rho1 + pf(ix1:ix2,jy1:jy2,kz1:kz2)*(1./rho2 - 1./rho1)

!           xnrom(ix1:ix2,jy1:jy2,kz1:kz2) * (rhof(ix1+1:ix2+1,jy1:jy2,kz1:kz2) -   &
!            rhof(ix1-1:ix2-1,jy1:jy2,kz1:kz2) )/dx + &
!           ynorm(ix1:ix2,jy1:jy2,kz1:kz2) * (rhof(ix1:ix2,jy1+1:jy2+1,kz1:kz2) -   &
!            rhof(ix1:ix2,jy1-1:jy2-1,kz1:kz2) )/dy


           aixr = 0.
           aixl = 0.
           aiyr = 0.
           aiyl = 0.
        testval = 0.


      do i=ix1,ix2
      do j=jy1,jy2
      do k=kz1,kz2

      !rhoxr = (smrh(i+1,j,k)+smrh(i,j,k))/2.d0
      !rhoxl = (smrh(i-1,j,k)+smrh(i,j,k))/2.d0
      !rhoyr = (smrh(i,j+1,k)+smrh(i,j,k))/2.d0
      !rhoyl = (smrh(i,j-1,k)+smrh(i,j,k))/2.d0
      
      rhoxr = 2.d0/(smrh(i+1,j,k)+smrh(i,j,k))
      rhoxl = 2.d0/(smrh(i-1,j,k)+smrh(i,j,k))
      rhoyr = 2.d0/(smrh(i,j+1,k)+smrh(i,j,k))
      rhoyl = 2.d0/(smrh(i,j-1,k)+smrh(i,j,k))

      !rhoxr = ((rho1x(i+1,j,k)+rho2x(i+1,j,k))/rho2)
      !rhoxl = ((rho1x(i,j,k)+rho2x(i,j,k))/rho2)
      !rhoyr = ((rho1y(i,j+1,k)+rho2y(i,j+1,k))/rho2)
      !rhoyl = ((rho1y(i,j,k)+rho2y(i,j,k))/rho2)

      aixr = (mdot(i,j,k)*xnorm(i,j,k)) * (rhoxr - 1./smrh(i,j,k))
      aixl = (mdot(i,j,k)*xnorm(i,j,k)) * (rhoxl - 1./smrh(i,j,k))
      aiyr = (mdot(i,j,k)*ynorm(i,j,k)) * (rhoyr - 1./smrh(i,j,k))
      aiyl = (mdot(i,j,k)*ynorm(i,j,k)) * (rhoyl - 1./smrh(i,j,k))

      !aixr = (mdot(i,j,k)+mdot(i+1,j,k)) * (xnorm(i,j,k)+xnorm(i+1,j,k)) * 0.25d0 * (rhoxr - 1./smrh(i,j,k))
      !aixl = (mdot(i,j,k)+mdot(i-1,j,k)) * (xnorm(i,j,k)+xnorm(i-1,j,k)) * 0.25d0 * (rhoxl - 1./smrh(i,j,k))
      !aiyr = (mdot(i,j,k)+mdot(i,j+1,k)) * (ynorm(i,j,k)+ynorm(i,j+1,k)) * 0.25d0 * (rhoyr - 1./smrh(i,j,k))
      !aiyl = (mdot(i,j,k)+mdot(i,j-1,k)) * (ynorm(i,j,k)+ynorm(i,j-1,k)) * 0.25d0 * (rhoyl - 1./smrh(i,j,k))

      !print *, "aixr",aixr-aixl
!      divv(i,j,k) =           & 
!         ( (uni(i+1,j,k) ) -   &
!           (uni(i,j,k)   ) )/dx +  &       
!         ( (vni(i,j+1,k) ) -   &
!           (vni(i,j,k)   ) )/dy

      testval(i,j,k) = (aixr-aixl)/dx + (aiyr-aiyl)/dy

      divv(i,j,k) =           & 
         ( (uni(i+1,j,k) + aixr) -   &
           (uni(i,j,k)   + aixl) )/dx +  &       
         ( (vni(i,j+1,k) + aiyr) -   &
           (vni(i,j,k)   + aiyl) )/dy

     ! divv(i,j,k) =           &
     !    (uni(i+1,j,k) - uni(i,j,k))/dx +  &
     !    (vni(i,j+1,k) - vni(i,j,k))/dy +  &
     !    (aixr         - aixl      )/dx +  &
     !    (aiyr         - aiyl      )/dy





      end do
      end do
      end do

      !rhoc = (rho_hxl(i,j)+rho_hxl(i+1,j)+rho_hyl(i,j)+rho_hyl(i,j+1))/4.d0
!      rhoxr = 2./(smrh(ix1+1,jy1,kz1)+smrh(ix1,jy1,kz1))
!      rhoxl = 2./(smrh(ix1-1,jy1,kz1)+smrh(ix1,jy1,kz1))
!      rhoyr = 2./(smrh(ix1,jy1+1,kz1)+smrh(ix1,jy1,kz1))
!      rhoyl = 2./(smrh(ix1,jy1-1,kz1)+smrh(ix1,jy1,kz1))
!      aixr = (mdot(ix1,jy1,kz1)*xnorm(ix1,jy1,kz1)) * (rhoxr - 1./smrh(ix1,jy1,kz1))
!      aixl = (mdot(ix1,jy1,kz1)*xnorm(ix1,jy1,kz1)) * (rhoxl - 1./smrh(ix1,jy1,kz1)) 
!      aiyr = (mdot(ix1,jy1,kz1)*ynorm(ix1,jy1,kz1)) * (rhoyr - 1./smrh(ix1,jy1,kz1))
!      aiyl = (mdot(ix1,jy1,kz1)*ynorm(ix1,jy1,kz1)) * (rhoyl - 1./smrh(ix1,jy1,kz1))


!      divv(ix1:ix2,jy1:jy2,kz1:kz2) =           & 
!         ( (uni(ix1+1:ix2+1,jy1:jy2,kz1:kz2) + aixr) -   &
!           (uni(ix1:ix2,jy1:jy2,kz1:kz2)     + aixl) )/dx +  &       
!         ( (vni(ix1:ix2,jy1+1:jy2+1,kz1:kz2) + aiyr) -   &
!           (vni(ix1:ix2,jy1:jy2,kz1:kz2)     + aiyl) )/dy
 
!      divv(ix1:ix2,jy1:jy2,kz1:kz2) =           & 
!         ( uni(ix1+1:ix2+1,jy1:jy2,kz1:kz2) -   &
!           uni(ix1:ix2,jy1:jy2,kz1:kz2) )/dx +  &       
!         ( vni(ix1:ix2,jy1+1:jy2+1,kz1:kz2) -   &
!           vni(ix1:ix2,jy1:jy2,kz1:kz2) )/dy + &
!          ( xnorm(ix1:ix2,jy1:jy2,kz1:kz2) * ( (rho2 + (rho1 - rho2) * smhv(ix1+1:ix2+1,jy1:jy2,kz1:kz2)) -   &
!            (rho2 + (rho1 - rho2) * smhv(ix1-1:ix2-1,jy1:jy2,kz1:kz2)) )/2./dx + &
!            ynorm(ix1:ix2,jy1:jy2,kz1:kz2) * ( (rho2 + (rho1 - rho2) * smhv(ix1:ix2,jy1+1:jy2+1,kz1:kz2))  -   &
!            (rho2 + (rho1 - rho2) * smhv(ix1:ix2,jy1-1:jy2-1,kz1:kz2)) )/2./dy ) 

! added by mslee
!     divv(ix1:ix2,jy1:jy2,kz1:kz2) =           & 
!         ( uni(ix1+1:ix2+1,jy1:jy2,kz1:kz2) -   &
!           uni(ix1:ix2,jy1:jy2,kz1:kz2) )/dx +  &       
!         ( vni(ix1:ix2,jy1+1:jy2+1,kz1:kz2) -   &
!           vni(ix1:ix2,jy1:jy2,kz1:kz2) )/dy


#if NDIM == 3
      divv(ix1:ix2,jy1:jy2,kz1:kz2) =           & 
           divv(ix1:ix2,jy1:jy2,kz1:kz2)     +  &
         ( wni(ix1:ix2,jy1:jy2,kz1+1:kz2+1) -   &
           wni(ix1:ix2,jy1:jy2,kz1:kz2) )/dz 
#endif

       ! deallocate(xnorm,ynorm,rhof)


END SUBROUTINE ins_divergence_VD_ML



SUBROUTINE ins_corrector_VD(uni,vni,wni,sigx,sigy,sigz,p,ix1,ix2,jy1,jy2,kz1,kz2, &
        dt,dx,dy,dz,alfa,rho1x,rho2x,rho1y,rho2y,rho1z,rho2z)

      ! This routine computes the corrected divergence-free velocities.
    
      use Grid_data,        ONLY : gr_meshMe

      implicit none

#include "Flash.h"

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dt,dx,dy,dz,alfa
      REAL, DIMENSION(:,:,:), INTENT(IN) :: p,rho1x,rho2x,rho1y,rho2y,rho1z,rho2z
      REAL, DIMENSION(:,:,:), INTENT(IN OUT) :: uni,vni,wni,sigx,sigy,sigz

      REAL :: coef,Mdens
      INTEGER :: i,j,k

      coef = dt*alfa

      !- kpd - Doesn't loop through boundary X-mom locations.
      !        Those were taken care of in ab2rk3
      do k=kz1,kz2
         do j=jy1,jy2
            do i=ix1+1,ix2

               !- kpd - inverse of the mixture density
               !--------------------------------------
               Mdens = (rho1x(i,j,k) + rho2x(i,j,k))
               !Mdens = 1.0 !(rho1x(i,j,k) + rho2x(i,j,k))
 
               !-------------------------------------------------
               !- kpd - Correct the final velocity...
               !           u(n+1) = u(*) - dt/rho*[ dP/dx-sig*K ]
               !        
               !    *** Note: sigx already contains 1/rho(mix)
               !-------------------------------------------------
               uni(i,j,k) =                                & 
                           uni(i,j,k) -                    &
                           coef*( Mdens*( p(i,j,k) -       &
                                          p(i-1,j,k) )/dx  &
                                  -sigx(i,j,k) )

               !if (gr_meshMe .eq. 0) print*,"KPD CorrX1",i,j,k,Mdens,sigx(i,j,k),uni(i,j,k)

            enddo
         enddo
      enddo

!      uni(ix1+1:ix2,jy1:jy2,kz1:kz2) =             & 
!         uni(ix1+1:ix2,jy1:jy2,kz1:kz2) -          &
!         coef*( ( p(ix1+1:ix2,jy1:jy2,kz1:kz2) -     &
!                  p(ix1:ix2-1,jy1:jy2,kz1:kz2) )/dx  &
!               -sigx(ix1+1:ix2,jy1:jy2,kz1:kz2))


      !- kpd - Doesn't loop through boundary Y-mom locations.
      !        Those were taken care of in ab2rk3
      do k=kz1,kz2
         do j=jy1+1,jy2
            do i=ix1,ix2

               !- kpd - inverse of the mixture density
               !--------------------------------------
               Mdens = (rho1y(i,j,k) + rho2y(i,j,k))
               !Mdens = 1.0 !(rho1y(i,j,k) + rho2y(i,j,k))

               !print*,"KPD CorrY",i,j,k,Mdens,sigy(i,j,k)

               !-------------------------------------------------
               !- kpd - Correct the final velocity...
               !           u(n+1) = u(*) - dt/rho*[ dP/dx-sig*K ]
               !        
               !        Note: sigy already contains 1/rho(mix)
               !-------------------------------------------------
!print*,"CORRECT",i,j,vni(i,j,k),Mdens,(p(i,j,k)-p(i,j-1,k) )/dy,coef*Mdens*(p(i,j,k)-p(i,j-1,k) )/dy
               vni(i,j,k) =                                &           
                           vni(i,j,k) -                    &
                           coef*( Mdens*( p(i,j,k) -       &
                                          p(i,j-1,k) )/dy  &
                                  -sigy(i,j,k) )

               !if (gr_meshMe .eq. 0) print*,"KPD CorrY1",i,j,k,Mdens,sigy(i,j,k),vni(i,j,k)
            enddo
         enddo
      enddo

!      vni(ix1:ix2,jy1+1:jy2,kz1:kz2) =             &
!         vni(ix1:ix2,jy1+1:jy2,kz1:kz2) -          &
!         coef*( ( p(ix1:ix2,jy1+1:jy2,kz1:kz2) -     &
!                  p(ix1:ix2,jy1:jy2-1,kz1:kz2) )/dy  &
!               -sigy(ix1:ix2,jy1+1:jy2,kz1:kz2))

#if NDIM == 3
!      wni(ix1:ix2,jy1:jy2,kz1+1:kz2) =             &
!         wni(ix1:ix2,jy1:jy2,kz1+1:kz2) -          &
!         coef*( ( p(ix1:ix2,jy1:jy2,kz1+1:kz2) -     &
!                  p(ix1:ix2,jy1:jy2,kz1:kz2-1) )/dz  &
!               -sigz(ix1:ix2,jy1:jy2,kz1+1:kz2))

      do k=kz1+1,kz2
         do j=jy1,jy2
            do i=ix1,ix2

               !- kpd - inverse of the mixture density
               Mdens = (rho1z(i,j,k) + rho2z(i,j,k))
               !Mdens = 1.0 !(rho1z(i,j,k) + rho2z(i,j,k))

               !-------------------------------------------------
               !- kpd - Correct the final velocity...
               !           u(n+1) = u(*) - dt/rho*[ dP/dx-sig*K ]
               !        
               !        Note: sigz already contains 1/rho(mix)
               !-------------------------------------------------
               wni(i,j,k) =                                &           
                           wni(i,j,k) -                    &
                           coef*( Mdens*( p(i,j,k) -       &
                                          p(i,j,k-1) )/dz  &
                           -sigz(i,j,k))
            enddo
         enddo
      enddo

#endif

END SUBROUTINE ins_corrector_VD

