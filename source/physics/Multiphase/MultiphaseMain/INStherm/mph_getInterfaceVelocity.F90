    SUBROUTINE mph_getInterfaceVelocity(uni,vni,ru1,ix1,ix2,jy1,jy2,dx,dy, &
                           visc,rho1x,rho2x,rho1y,rho2y,gravX,gravY, &
                           mdot,smrh,xnorm,ynorm,uint,vint)

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
      INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2
      REAL, INTENT(IN):: ru1, dx, dy, gravX,gravY
      REAL, DIMENSION(:,:,:), INTENT(IN):: uni, vni, visc, rho1x, rho2x, rho1y, rho2y
      REAL, DIMENSION(:,:,:), INTENT(IN) :: xnorm,ynorm,mdot,smrh
      REAL, DIMENSION(:,:,:), INTENT(OUT):: uint, vint

      INTEGER:: i, j
      REAL:: dx1, dy1, Mdens
      INTEGER, parameter :: kz1 = 1


      !ML - ghost fluid for u and v
      INTEGER:: ig, jg
      REAL:: rhox,rhoy,mdotx,mdoty,normx,normy,rhoxc,rhoyc 


      !--- Interface Velocity Method 1 ---!

      !++++++++++  U-COMPONENT  ++++++++++

       !uint(ix1:ix2+1,jy1:jy2,kz1) =  uni(ix1:ix2+1,jy1:jy2,kz1) + &
       !                               (mdot(ix1-1:ix2,jy1:jy2,kz1) + mdot(ix1:ix2+1,jy1:jy2,kz1))/2.0d0 * &
       !                               (xnorm(ix1-1:ix2,jy1:jy2,kz1) + xnorm(ix1:ix2+1,jy1:jy2,kz1))/2.0d0 * &
       !                               (rho1x(ix1:ix2+1,jy1:jy2,kz1) + rho2x(ix1:ix2+1,jy1:jy2,kz1))

      !++++++++++  V-COMPONENT  ++++++++++

       !vint(ix1:ix2,jy1:jy2+1,kz1) =  vni(ix1:ix2,jy1:jy2+1,kz1) + &
       !                               (mdot(ix1:ix2,jy1-1:jy2,kz1) + mdot(ix1:ix2,jy1:jy2+1,kz1))/2.0d0 * &
       !                               (ynorm(ix1:ix2,jy1-1:jy2,kz1) + ynorm(ix1:ix2,jy1:jy2+1,kz1))/2.0d0 * &
       !                               (rho1y(ix1:ix2,jy1:jy2+1,kz1) + rho2y(ix1:ix2,jy1:jy2+1,kz1))


      !--- Interface Velocity Method 2 ---!

      !++++++++++  U-COMPONENT  ++++++++++

       uint(ix1:ix2+1,jy1:jy2,kz1) =  uni(ix1:ix2+1,jy1:jy2,kz1) + &
                                      (mdot(ix1-1:ix2,jy1:jy2,kz1) + mdot(ix1:ix2+1,jy1:jy2,kz1))/2.0d0 * &
                                      (xnorm(ix1-1:ix2,jy1:jy2,kz1) + xnorm(ix1:ix2+1,jy1:jy2,kz1))/2.0d0 * &
                                      (smrh(ix1-1:ix2,jy1:jy2,kz1) + smrh(ix1:ix2+1,jy1:jy2,kz1))/2.0d0

      !++++++++++  V-COMPONENT  ++++++++++

       vint(ix1:ix2,jy1:jy2+1,kz1) =  vni(ix1:ix2,jy1:jy2+1,kz1) + &
                                      (mdot(ix1:ix2,jy1-1:jy2,kz1) + mdot(ix1:ix2,jy1:jy2+1,kz1))/2.0d0 * &
                                      (ynorm(ix1:ix2,jy1-1:jy2,kz1) + ynorm(ix1:ix2,jy1:jy2+1,kz1))/2.0d0 * &
                                      (smrh(ix1:ix2,jy1-1:jy2,kz1) + smrh(ix1:ix2,jy1:jy2+1,kz1))/2.0d0


       END SUBROUTINE mph_getInterfaceVelocity

