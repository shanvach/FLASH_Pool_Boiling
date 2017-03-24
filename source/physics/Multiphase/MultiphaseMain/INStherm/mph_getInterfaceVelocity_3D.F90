    SUBROUTINE mph_getInterfaceVelocity_3D(uni,vni,wni,ru1,ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz,&
                           visc,rho1x,rho2x,rho1y,rho2y,rho1z,rho2z,gravX,gravY,gravZ,&
                           mdot,smrh,xnorm,ynorm,znorm,uint,vint,wint)

  !-------------------------------------------------------------
  ! Computes the interface velocity to be used in level 
  ! advection from the fluid velocity.
  ! 
  !  Input: uni,vni      = velocity at timestep n+1
  !         mdot         = mass flux
  !         smrh         = smeared density
  !         xnorm,ynorm  = interface normal
  !         jy1,jy2      = starting and ending y indices
  !         dx,dy        = grid spacing in x and y directions
  !  
  ! Output: uint,vint    = interface velocity at timestep n+1
  !-------------------------------------------------------------

      use Driver_interface, ONLY : Driver_abortFlash

      use Driver_data,      ONLY : dr_nstep

      use RuntimeParameters_interface, ONLY : RuntimeParameters_get

      use IncompNS_data, ONLY : ins_iConvU

      use Multiphase_data, ONLY : mph_rho2

#include "Flash.h"
!#include "constants.h"

      implicit none
      INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2, kz1, kz2
      REAL, INTENT(IN):: ru1, dx, dy, dz, gravX, gravY, gravZ
      REAL, DIMENSION(:,:,:), INTENT(IN):: uni, vni, wni, visc, rho1x, rho2x, rho1y, rho2y, rho1z, rho2z
      REAL, DIMENSION(:,:,:), INTENT(IN) :: xnorm,ynorm,znorm,mdot,smrh
      REAL, DIMENSION(:,:,:), INTENT(OUT):: uint, vint, wint

      INTEGER:: i, j, k
      REAL:: dx1, dy1, Mdens
      !ML - ghost fluid for u and v
      INTEGER:: ig, jg
      REAL:: rhox,rhoy,mdotx,mdoty,normx,normy 
      REAL:: rhoz,mdotz,normz


      !--- Interface Velocity Method 1 ---!

      !____________________________U component________________________________________!

      uint(ix1:ix2+1,jy1:jy2,kz1:kz2) = uni(ix1:ix2+1,jy1:jy2,kz1:kz2) + &
                                      (mdot(ix1-1:ix2,jy1:jy2,kz1:kz2) + mdot(ix1:ix2+1,jy1:jy2,kz1:kz2))/2.0d0 * &
                                      (xnorm(ix1-1:ix2,jy1:jy2,kz1:kz2) + xnorm(ix1:ix2+1,jy1:jy2,kz1:kz2))/2.0d0 * &
                                      (rho1x(ix1:ix2+1,jy1:jy2,kz1:kz2) + rho2x(ix1:ix2+1,jy1:jy2,kz1:kz2))


      !____________________________V component________________________________________!

      vint(ix1:ix2,jy1:jy2+1,kz1:kz2) = vni(ix1:ix2,jy1:jy2+1,kz1:kz2) + &
                                      (mdot(ix1:ix2,jy1-1:jy2,kz1:kz2) + mdot(ix1:ix2,jy1:jy2+1,kz1:kz2))/2.0d0 * &
                                      (ynorm(ix1:ix2,jy1-1:jy2,kz1:kz2) + ynorm(ix1:ix2,jy1:jy2+1,kz1:kz2))/2.0d0 * &
                                      (rho1y(ix1:ix2,jy1:jy2+1,kz1:kz2) + rho2y(ix1:ix2,jy1:jy2+1,kz1:kz2))

      !____________________________W component________________________________________!

      wint(ix1:ix2,jy1:jy2,kz1:kz2+1) = wni(ix1:ix2,jy1:jy2,kz1:kz2+1) + &
                                      (mdot(ix1:ix2,jy1:jy2,kz1-1:kz2) + mdot(ix1:ix2,jy1:jy2,kz1:kz2+1))/2.0d0 * &
                                      (znorm(ix1:ix2,jy1:jy2,kz1-1:kz2) + znorm(ix1:ix2,jy1:jy2,kz1:kz2+1))/2.0d0 * &
                                      (rho1z(ix1:ix2,jy1:jy2,kz1:kz2+1) + rho2z(ix1:ix2,jy1:jy2,kz1:kz2+1))




      !--- Interface Velocity Method 2 ---!

      !____________________________U component________________________________________!

      !uint(ix1:ix2+1,jy1:jy2,kz1:kz2) = uni(ix1:ix2+1,jy1:jy2,kz1:kz2) + &
      !                                (mdot(ix1-1:ix2,jy1:jy2,kz1:kz2) + mdot(ix1:ix2+1,jy1:jy2,kz1:kz2))/2.0d0 * &
      !                                (xnorm(ix1-1:ix2,jy1:jy2,kz1:kz2) + xnorm(ix1:ix2+1,jy1:jy2,kz1:kz2))/2.0d0 * &
      !                                (smrh(ix1-1:ix2,jy1:jy2,kz1:kz2) + smrh(ix1:ix2+1,jy1:jy2,kz1:kz2))/2.0d0


      !____________________________V component________________________________________!

      !vint(ix1:ix2,jy1:jy2+1,kz1:kz2) = vni(ix1:ix2,jy1:jy2+1,kz1:kz2) + &
      !                                (mdot(ix1:ix2,jy1-1:jy2,kz1:kz2) + mdot(ix1:ix2,jy1:jy2+1,kz1:kz2))/2.0d0 * &
      !                                (ynorm(ix1:ix2,jy1-1:jy2,kz1:kz2) + ynorm(ix1:ix2,jy1:jy2+1,kz1:kz2))/2.0d0 * &
      !                                (smrh(ix1:ix2,jy1-1:jy2,kz1:kz2) + smrh(ix1:ix2,jy1:jy2+1,kz1:kz2))/2.0d0

      !____________________________W component________________________________________!

      !wint(ix1:ix2,jy1:jy2,kz1:kz2+1) = wni(ix1:ix2,jy1:jy2,kz1:kz2+1) + &
      !                                (mdot(ix1:ix2,jy1:jy2,kz1-1:kz2) + mdot(ix1:ix2,jy1:jy2,kz1:kz2+1))/2.0d0 * &
      !                                (znorm(ix1:ix2,jy1:jy2,kz1-1:kz2) + znorm(ix1:ix2,jy1:jy2,kz1:kz2+1))/2.0d0 * &
      !                                (smrh(ix1:ix2,jy1:jy2,kz1-1:kz2) + smrh(ix1:ix2,jy1:jy2,kz1:kz2+1))/2.0d0


END SUBROUTINE mph_getInterfaceVelocity_3D

