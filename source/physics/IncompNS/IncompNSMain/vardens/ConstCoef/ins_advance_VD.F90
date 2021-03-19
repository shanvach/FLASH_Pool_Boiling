!########################################################################

SUBROUTINE ins_predictor_VD(uni,vni,wni,unew,vnew,wnew,uold,vold,&
        wold,p,dt,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,gama,rhoa,alfa)

      ! This routine computes the intermediate velocities based on
      ! the explicit second-order Adams-Bashforth scheme (gamma=1.5,
      ! rhoa=-0.5,alfa=1), or a third order Runge-Kutta method.

      use IncompNS_data, ONLY : ins_prescoeff

      implicit none

#include "Flash.h"

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dt
      REAL, DIMENSION(:), INTENT(IN) :: dx,dy,dz
      REAL, DIMENSION(:,:,:), INTENT(IN) :: unew,vnew,wnew,&
                                            uold,vold,wold,&
                                            p
      REAL, DIMENSION(:,:,:), INTENT(IN OUT) :: uni,vni,wni

      REAL :: gama,rhoa,alfa
      INTEGER :: i,j,k

      !---------------------------------------------------------------------
      !---------------------------------------------------------------------
      !- kpd - When pressure correction scheme is not used ins_prescoeff = 0
      !           and dP/dx is not used in u* predictor step
      !        uni and vni are u* and v*
      !---------------------------------------------------------------------
      !---------------------------------------------------------------------

      do i=ix1,ix2+1

        uni(i,jy1:jy2,kz1:kz2) = uni(i,jy1:jy2,kz1:kz2) + &
           dt*(gama*unew(i,jy1:jy2,kz1:kz2) +                     &
               rhoa*uold(i,jy1:jy2,kz1:kz2) -                     &
               ins_prescoeff*alfa*( p(i,jy1:jy2,kz1:kz2) -        &
                                    p(i-1,jy1:jy2,kz1:kz2))*dx(i))
      enddo

      do j=jy1,jy2+1

        vni(ix1:ix2,j,kz1:kz2) = vni(ix1:ix2,j,kz1:kz2) + &
           dt*(gama*vnew(ix1:ix2,j,kz1:kz2) +                     &
               rhoa*vold(ix1:ix2,j,kz1:kz2) -                     &
               ins_prescoeff*alfa*( p(ix1:ix2,j,kz1:kz2) -        &
                                    p(ix1:ix2,j-1,kz1:kz2) )*dy(j))
      enddo

#if NDIM == 3

      do k=kz1,kz2+1

        wni(ix1:ix2,jy1:jy2,k) = wni(ix1:ix2,jy1:jy2,k) + &
           dt*(gama*wnew(ix1:ix2,jy1:jy2,k) +                     &
               rhoa*wold(ix1:ix2,jy1:jy2,k) -                     &
               ins_prescoeff*alfa*( p(ix1:ix2,jy1:jy2,k) -        &
                                    p(ix1:ix2,jy1:jy2,k-1) )*dz(k))
      enddo

#endif


END SUBROUTINE ins_predictor_VD

!########################################################################

SUBROUTINE ins_predictor_VD_IB(uni,vni,wni,unew,vnew,wnew,uold,vold,&
        wold,p,dt,ix1,ix2,jy1,jy2,kz1,kz2,gama,rhoa,alfa,poldx,poldy,poldz)

      ! This routine computes the intermediate velocities based on
      ! the explicit second-order Adams-Bashforth scheme (gamma=1.5,
      ! rhoa=-0.5,alfa=1), or a third order Runge-Kutta method.

      use IncompNS_data, ONLY : ins_prescoeff

      implicit none

#include "Flash.h"

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dt
      REAL, DIMENSION(:,:,:), INTENT(IN) :: unew,vnew,wnew,&
                                            uold,vold,wold,&
                                            p,poldx,poldy,poldz
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
                 poldx(ix1:ix2+1,jy1:jy2,kz1:kz2))

      vni(ix1:ix2,jy1:jy2+1,kz1:kz2) = vni(ix1:ix2,jy1:jy2+1,kz1:kz2) + &
         dt*(gama*vnew(ix1:ix2,jy1:jy2+1,kz1:kz2) +                     &
             rhoa*vold(ix1:ix2,jy1:jy2+1,kz1:kz2) -                     &
                 poldy(ix1:ix2,jy1:jy2+1,kz1:kz2))

#if NDIM == 3
      wni(ix1:ix2,jy1:jy2,kz1:kz2+1) = wni(ix1:ix2,jy1:jy2,kz1:kz2+1) + &
         dt*(gama*wnew(ix1:ix2,jy1:jy2,kz1:kz2+1) +                     &
             rhoa*wold(ix1:ix2,jy1:jy2,kz1:kz2+1) -                     &
                 poldz(ix1:ix2,jy1:jy2,kz1:kz2+1))
#endif


END SUBROUTINE ins_predictor_VD_IB

!########################################################################

SUBROUTINE ins_divergence_VD(uni,vni,wni,ix1,ix2,jy1,jy2,kz1,kz2,&
         dx,dy,dz,divv)

      ! This routine computes the divergence of the velocity field.

      implicit none

#include "Flash.h"

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, DIMENSION(:), INTENT(IN) :: dx,dy,dz
      REAL, DIMENSION(:,:,:), INTENT(IN) :: uni,vni,wni
      REAL, DIMENSION(:,:,:), INTENT(OUT) :: divv
      
      INTEGER :: i,j,k

      do i=ix1,ix2
         do j=jy1,jy2

            divv(i,j,kz1:kz2) =           & 
                ( uni(i+1,j  ,kz1:kz2) -   &
                  uni(i  ,j  ,kz1:kz2) )*dx(i) +  &       
                ( vni(i  ,j+1,kz1:kz2) -   &
                  vni(i  ,j  ,kz1:kz2) )*dy(j)
         enddo
      enddo

#if NDIM == 3
     do k=kz1,kz2

      divv(ix1:ix2,jy1:jy2,k) =           & 
           divv(ix1:ix2,jy1:jy2,k)     +  &
         ( wni(ix1:ix2,jy1:jy2,k+1) -   &
           wni(ix1:ix2,jy1:jy2,k  ) )*dz(k) 
     enddo
#endif

END SUBROUTINE ins_divergence_VD

!########################################################################

SUBROUTINE ins_firstCorrector_VD_IB(uni,vni,wni,ix1,ix2,jy1,jy2,kz1,kz2,dt,poldx,poldy,poldz)

      ! This routine computes the corrected divergence-free velocities.
    
      use Grid_data,        ONLY : gr_meshMe

      implicit none

#include "Flash.h"

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dt
      REAL, DIMENSION(:,:,:), INTENT(IN OUT) :: uni,vni,wni
      REAL, DIMENSION(:,:,:), INTENT(IN) :: poldx,poldy,poldz

      INTEGER :: i,j,k
      REAL :: coef

      coef = dt

      !- kpd - Doesn't loop through boundary X-mom locations.
      !        Those were taken care of in ab2rk3
      do k=kz1,kz2
         do j=jy1,jy2
            do i=ix1,ix2+1
               uni(i,j,k) = uni(i,j,k) + coef*poldx(i,j,k)
            enddo
         enddo
      enddo

      do k=kz1,kz2
         do j=jy1,jy2+1
            do i=ix1,ix2
               vni(i,j,k) = vni(i,j,k) + coef*poldy(i,j,k)
            enddo
         enddo
      enddo

#if NDIM == 3
      do k=kz1,kz2+1
         do j=jy1,jy2
            do i=ix1,ix2
               wni(i,j,k) = wni(i,j,k) + coef*poldz(i,j,k)
            enddo
         enddo
      enddo
#endif

END SUBROUTINE ins_firstCorrector_VD_IB

