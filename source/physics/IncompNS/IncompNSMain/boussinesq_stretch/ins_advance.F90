SUBROUTINE ins_predictor(uni,vni,wni,unew,vnew,wnew,uold,vold, &
                         wold,p,dt,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,gama,rhoa,alfa)

      ! This routine computes the intermediate velocities based on
      ! the explicit second-order Adams-Bashforth scheme (gamma=1.5,
      ! rhoa=-0.5,alfa=1), or a third order Runge-Kutta method.

      use IncompNS_data, ONLY : ins_prescoeff,ins_dpdx,ins_dpdy,ins_dpdz, &
                                ins_gravX,ins_gravY,ins_gravZ
      
      use Simulation_data, ONLY : sim_xMin, sim_yMin, sim_zMin, &
                                  sim_xMax, sim_yMax, sim_zMax

      implicit none

#include "Flash.h"
#include "constants.h"

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dt, gama,rhoa,alfa
      REAL, DIMENSION(:), INTENT(IN) :: dx,dy,dz
      REAL, DIMENSION(:,:,:), INTENT(IN) :: unew,vnew,wnew,uold,vold,wold,p
      REAL, DIMENSION(:,:,:), INTENT(IN OUT) :: uni,vni,wni

      INTEGER :: i,j,k
      REAL    :: dpx, dpy, dpz

      ! calculate delta press to correct for stretched grid
      dpx = ins_dpdx * (sim_xMax - sim_xMin)
      dpy = ins_dpdy * (sim_yMax - sim_yMin)
      dpz = ins_dpdz * (sim_zMax - sim_zMin)

      do k = kz1,kz2
        do j = jy1,jy2
          uni(ix1:ix2+1,j,k) = uni(ix1:ix2+1,j,k) + dt*(                             &
            gama*unew(ix1:ix2+1,j,k) + rhoa*uold(ix1:ix2+1,j,k) -                    &
            ins_prescoeff*alfa*dx(ix1:ix2+1)*(p(ix1:ix2+1,j,k) - p(ix1-1:ix2,j,k)) - &
            !alfa*dpx*dx(ix1:ix2+1) + alfa*ins_gravX )
            alfa*ins_dpdx + alfa*ins_gravX )
        end do
      end do

      do k = kz1,kz2
        do i = ix1,ix2
          vni(i,jy1:jy2+1,k) = vni(i,jy1:jy2+1,k) + dt*(                             &
            gama*vnew(i,jy1:jy2+1,k) + rhoa*vold(i,jy1:jy2+1,k) -                    &
            ins_prescoeff*alfa*dy(jy1:jy2+1)*(p(i,jy1:jy2+1,k) - p(i,jy1-1:jy2,k)) - &
            !alfa*dpy*dy(jy1:jy2+1) + alfa*ins_gravY )
            alfa*ins_dpdy + alfa*ins_gravY )
        end do
      end do    

#if NDIM == MDIM
      do j = jy1,jy2
        do i = ix1,ix2
          wni(i,j,kz1:kz2+1) = wni(i,j,kz1:kz2+1) + dt*(                             &
            gama*wnew(i,j,kz1:kz2+1) + rhoa*wold(i,j,kz1:kz2+1) -                    &
            ins_prescoeff*alfa*dz(kz1:kz2+1)*(p(i,j,kz1:kz2+1) - p(i,j,kz1-1:kz2)) - &
            !alfa*dpz*dz(kz1:kz2+1) + alfa*ins_gravZ )
            alfa*ins_dpdz + alfa*ins_gravZ )
        end do
      end do
#endif

END SUBROUTINE ins_predictor

!########################################################################

SUBROUTINE ins_corrector(uni,vni,wni,p,ix1,ix2,jy1,jy2,kz1,kz2, &
        dt,dx,dy,dz,alfa)

      ! This routine computes the corrected divergence-free velocities.
    
      implicit none

#include "Flash.h"
#include "constants.h"

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, INTENT(IN) :: dt,alfa
      REAL, DIMENSION(:), INTENT(IN) :: dx,dy,dz
      REAL, DIMENSION(:,:,:), INTENT(IN) :: p
      REAL, DIMENSION(:,:,:), INTENT(IN OUT) :: uni,vni,wni

      INTEGER :: i,j,k

      do k = kz1,kz2
        do j = jy1,jy2
          uni(ix1:ix2+1,j,k) = uni(ix1:ix2+1,j,k) - &
            dt*alfa*dx(ix1:ix2+1)*( p(ix1:ix2+1,j,k) - p(ix1-1:ix2,j,k) )
        end do
      end do

      do k = kz1,kz2
        do i = ix1,ix2
          vni(i,jy1:jy2+1,k) = vni(i,jy1:jy2+1,k) - &
            dt*alfa*dy(jy1:jy2+1)*( p(i,jy1:jy2+1,k) - p(i,jy1-1:jy2,k) )
        end do
      end do

#if NDIM == MDIM
      do j = jy1,jy2
        do i = ix1,ix2
          wni(i,j,kz1:kz2+1) = wni(i,j,kz1:kz2+1) - &
            dt*alfa*dz(kz1:kz2+1)*( p(i,j,kz1:kz2+1) - p(i,j,kz1-1:kz2) )
        end do
      end do
#endif

END SUBROUTINE ins_corrector

!########################################################################

SUBROUTINE ins_divergence(uni,vni,wni,ix1,ix2,jy1,jy2,kz1,kz2,&
                          dx,dy,dz,divv)

      ! This routine computes the divergence of the velocity field.

      implicit none

#include "Flash.h"
#include "constants.h"

      INTEGER, INTENT(IN) :: ix1,ix2,jy1,jy2,kz1,kz2
      REAL, DIMENSION(:), INTENT(IN) :: dx,dy,dz
      REAL, DIMENSION(:,:,:), INTENT(IN) :: uni,vni,wni
      REAL, DIMENSION(:,:,:), INTENT(OUT) :: divv

      INTEGER :: i, j, k

      do k = kz1,kz2
        do j = jy1,jy2
          do i = ix1,ix2
#if NDIM == MDIM
            divv(i,j,k) = dx(i) * (uni(i+1,j,k) - uni(i,j,k)) + &
                          dy(j) * (vni(i,j+1,k) - vni(i,j,k)) + &
                          dz(k) * (wni(i,j,k+1) - wni(i,j,k))
#else
            divv(i,j,k) = dx(i) * (uni(i+1,j,k) - uni(i,j,k)) + &
                          dy(j) * (vni(i,j+1,k) - vni(i,j,k))
#endif
          end do
        end do
      end do

END SUBROUTINE ins_divergence

!########################################################################
