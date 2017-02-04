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




      !++++++++++  U-COMPONENT  ++++++++++
        !do j = jy1,jy2
        !do i = ix1,ix2+1
        !     rhox  = 2.0d0/(smrh(i-1,j,kz1)+smrh(i,j,kz1))
             !rhox  = (smrh(i-1,j,kz1)+smrh(i,j,kz1))/2.d0
             !rhox =  (rho1x(i,j,kz1)+rho2x(i,j,kz1))/mph_rho2
        !     mdotx = (mdot(i-1,j,kz1)+mdot(i,j,kz1))/2.d0
        !     normx = (xnorm(i-1,j,kz1)+xnorm(i,j,kz1))/2.d0
             !normx = xnorm(i,j,kz1)
        !     uint(i,j,kz1)=uni(i,j,kz1) + mdotx*normx*rhox!/mph_rho2
        !end do
        !end do

        uint(ix1:ix2+1,jy1:jy2,kz1) = uni(ix1:ix2+1,jy1:jy2,kz1) + &
                                      (mdot(ix1-1:ix2,jy1:jy2,kz1) + mdot(ix1:ix2+1,jy1:jy2,kz1))/2.d0 * &
                                      (xnorm(ix1-1:ix2,jy1:jy1,kz1) + xnorm(ix1:ix2+1,jy1:jy2,kz1))/2.d0 * &
                                      (2.0d0/(smrh(ix1-1:ix2,jy1:jy2,kz1) + smrh(ix1:ix2+1,jy1:jy2,kz1)))

      !++++++++++  V-COMPONENT  ++++++++++
        !do j = jy1,jy2+1
        !do i = ix1,ix2
        !     rhoy  = 2.0d0/(smrh(i,j-1,kz1)+smrh(i,j,kz1))
             !rhoy  = (smrh(i,j-1,kz1)+smrh(i,j,kz1))/2.d0
             !rhoy  = (rho1y(i,j,kz1)+rho2y(i,j,kz1))/mph_rho2
        !     mdoty = (mdot(i,j-1,kz1)+mdot(i,j,kz1))/2.d0
        !     normy = (ynorm(i,j-1,kz1)+ynorm(i,j,kz1))/2.d0
             !normy = ynorm(i,j,kz1)
        !     vint(i,j,kz1)=vni(i,j,kz1) + mdoty*normy*rhoy!/mph_rho2
        !enddo
        !enddo

        vint(ix1:ix2+1,jy1:jy2,kz1) = vni(ix1:ix2+1,jy1:jy2,kz1) + &
                                      (mdot(ix1:ix2,jy1-1:jy2,kz1) + mdot(ix1:ix2,jy1:jy2+1,kz1))/2.d0 * &
                                      (ynorm(ix1:ix2,jy1-1:jy1,kz1) + ynorm(ix1:ix2,jy1:jy2+1,kz1))/2.d0 * &
                                      (2.0d0/(smrh(ix1:ix2,jy1-1:jy2,kz1) + smrh(ix1:ix2,jy1:jy2+1,kz1)))



       END SUBROUTINE mph_getInterfaceVelocity

