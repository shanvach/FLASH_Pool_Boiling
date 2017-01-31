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



      !++++++++++  U-COMPONENT  ++++++++++
        do k = kz1,kz2
        do j = jy1,jy2
        do i = ix1,ix2+1
             rhox  = 2.0d0/(smrh(i-1,j,kz1)+smrh(i,j,k))
             mdotx = (mdot(i-1,j,k)+mdot(i,j,k))/2.d0
             normx = (xnorm(i-1,j,k)+xnorm(i,j,k))/2.d0
             uint(i,j,k)=uni(i,j,k) + mdotx*normx*rhox
        end do
        end do
        end do

      !++++++++++  V-COMPONENT  ++++++++++
        do k = kz1,kz2
        do j = jy1,jy2+1
        do i = ix1,ix2
             rhoy  = 2.0d0/(smrh(i,j-1,k)+smrh(i,j,k))
             mdoty = (mdot(i,j-1,k)+mdot(i,j,k))/2.d0
             normy = (ynorm(i,j-1,k)+ynorm(i,j,k))/2.d0
             vint(i,j,k)=vni(i,j,k) + mdoty*normy*rhoy
        enddo
        enddo
        enddo

      !++++++++++  W-COMPONENT  ++++++++++
        do k = kz1,kz2+1
        do j = jy1,jy2
        do i = ix1,ix2
             rhoz  = 2.0d0/(smrh(i,j,k-1)+smrh(i,j,k))
             mdotz = (mdot(i,j,k-1)+mdot(i,j,k))/2.d0
             normz = (znorm(i,j,k-1)+znorm(i,j,k))/2.d0
             wint(i,j,k)=wni(i,j,k) + mdotz*normz*rhoz
        enddo
        enddo
        enddo

END SUBROUTINE mph_getInterfaceVelocity_3D

