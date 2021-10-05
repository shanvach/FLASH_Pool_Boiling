!=========================================================================
!=========================================================================
!=========================================================================


        subroutine mph_KPDcurvature3DAB(s,crv,dx,dy,dz, &
           ix1,ix2,jy1,jy2,kz1,kz2, &
           rho1x,rho2x,rho1y,rho2y,rho1z,rho2z,pf,rho1,rho2,visc,vis1,vis2)

        implicit none

#include "Flash.h"
#include "constants.h"

        !*****************************************************************

        !---------------------------------
        !- kpd - Data from routine call...
        !---------------------------------
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        real, intent(in) :: rho1, rho2, vis1, vis2
        real, dimension(:,:), intent(in) :: dx,dy,dz
        real, dimension(:,:,:), intent(inout):: s,crv, &
                                                rho1x,rho2x,rho1y, &
                                                rho2y,pf, &
                                                rho1z,rho2z, visc

        !--------------------------
        !- kpd - Local variables...
        !--------------------------
        integer :: i,j,k
        real :: rPhiXN,rPhiXE,rPhiXS,rPhiXW, &
                rPhiYN,rPhiYE,rPhiYS,rPhiYW, &
                rMagN,rMagE,rMagS,rMagW, &
                rPhiZN,rPhiZE,rPhiZS,rPhiZW, &
                rMagF,rMagB, &
                rPhiXF,rPhiXB,rPhiYF,rPhiYB,rPhiZF,rPhiZB

        real :: th,aa,xijl,xijr, &
                a1,a2,cri,xid,xij,xidl,xidr,yid,yidr,yidl,yij,yijl,yijr, &
                zijl,zijr,zid,zij,zidl,zidr
        real, parameter :: eps = 1E-13

        !*****************************************************************

        !---------------------------------------------
        !- kpd - Compute the curvature ---------------
        !---------------------------------------------

        crv = 0.
        do k = kz1,kz2
           do j = jy1,jy2
              do i = ix1,ix2    

              !-------------------------------------------------
              !- kpd - 3 phi gradients per face method
              !-------------------------------------------------

              !- kpd - Compute [d(phi)/dx] on all faces
              rPhiXE = dx(i,RIGHT_EDGE)*(s(i+1,j,k) - s(i,j,k)  )
              rPhiXW = dx(i,LEFT_EDGE)*(s(i,j,k)   - s(i-1,j,k))

              rPhiXN = (1./2)*(dx(i,  RIGHT_EDGE)*(s(i+1,j+1,k)-s(i  ,j+1,k))+&
                               dx(i-1,RIGHT_EDGE)*(s(i  ,j+1,k)-s(i-1,j+1,k))+&
                               dx(i,  RIGHT_EDGE)*(s(i+1,j  ,k)-s(i  ,j,k))+&
                               dx(i-1,RIGHT_EDGE)*(s(i  ,j  ,k)-s(i-1,j,k)))

              rPhiXS = (1./2)*(dx(i,  RIGHT_EDGE)*(s(i+1,j  ,k)-s(i  ,j,k))+&
                               dx(i-1,RIGHT_EDGE)*(s(i  ,j  ,k)-s(i-1,j,k))+&
                               dx(i,  RIGHT_EDGE)*(s(i+1,j-1,k)-s(i  ,j-1,k))+&
                               dx(i-1,RIGHT_EDGE)*(s(i  ,j-1,k)-s(i-1,j-1,k)))

              rPhiXF = (1./2)*(dx(i,  RIGHT_EDGE)*(s(i+1,j,k+1)-s(i  ,j,k+1))+&
                               dx(i-1,RIGHT_EDGE)*(s(i  ,j,k+1)-s(i-1,j,k+1))+&
                               dx(i,  RIGHT_EDGE)*(s(i+1,j,k  )-s(i  ,j,k  ))+&
                               dx(i-1,RIGHT_EDGE)*(s(i  ,j,k  )-s(i-1,j,k  )))

              rPhiXB = (1./2)*(dx(i,  RIGHT_EDGE)*(s(i+1,j,k  )-s(i  ,j,k  ))+&
                               dx(i-1,RIGHT_EDGE)*(s(i  ,j,k  )-s(i-1,j,k  ))+&
                               dx(i,  RIGHT_EDGE)*(s(i+1,j,k-1)-s(i  ,j,k-1))+&
                               dx(i-1,RIGHT_EDGE)*(s(i  ,j,k-1)-s(i-1,j,k-1)))

              !- kpd - Compute [d(phi)/dy] on all faces

              rPhiYE = (1./2)*(dy(j,  RIGHT_EDGE)*(s(i+1,j+1,k)-s(i+1,j  ,k))+&
                               dy(j-1,RIGHT_EDGE)*(s(i+1,j  ,k)-s(i+1,j-1,k))+&
                               dy(j,  RIGHT_EDGE)*(s(i  ,j+1,k)-s(i  ,j  ,k))+&
                               dy(j-1,RIGHT_EDGE)*(s(i  ,j  ,k)-s(i  ,j-1,k)))

              rPhiYW = (1./2)*(dy(j,  RIGHT_EDGE)*(s(i  ,j+1,k)-s(i  ,j  ,k))+&
                               dy(j-1,RIGHT_EDGE)*(s(i  ,j  ,k)-s(i  ,j-1,k))+&
                               dy(j,  RIGHT_EDGE)*(s(i-1,j+1,k)-s(i-1,j  ,k))+&
                               dy(j-1,RIGHT_EDGE)*(s(i-1,j  ,k)-s(i-1,j-1,k)))

              rPhiYN = dy(j,RIGHT_EDGE)*(s(i,j+1,k) - s(i,j,k)  )
              rPhiYS = dy(j,LEFT_EDGE)*(s(i,j,k)   - s(i,j-1,k))


              rPhiYF = (1./2)*(dy(j,  RIGHT_EDGE)*(s(i,j+1,k+1)-s(i,j  ,k+1))+&
                               dy(j-1,RIGHT_EDGE)*(s(i,j  ,k+1)-s(i,j-1,k+1))+&
                               dy(j,  RIGHT_EDGE)*(s(i,j+1,k  )-s(i,j  ,k  ))+&
                               dy(j-1,RIGHT_EDGE)*(s(i,j  ,k  )-s(i,j-1,k  )))

              rPhiYB = (1./2)*(dy(j,  RIGHT_EDGE)*(s(i,j+1,k  )-s(i,j  ,k  ))+&
                               dy(j-1,RIGHT_EDGE)*(s(i,j  ,k  )-s(i,j-1,k  ))+&
                               dy(j,  RIGHT_EDGE)*(s(i,j+1,k-1)-s(i,j  ,k-1))+&
                               dy(j-1,RIGHT_EDGE)*(s(i,j  ,k-1)-s(i,j-1,k-1)))


              !- kpd - Compute [d(phi)/dz] on all faces

              rPhiZE = (1./2)*(dz(k,  RIGHT_EDGE)*(s(i+1,j,k+1)-s(i+1,j,k  ))+&
                               dz(k-1,RIGHT_EDGE)*(s(i+1,j  ,k)-s(i+1,j,k-1))+&
                               dz(k,  RIGHT_EDGE)*(s(i  ,j,k+1)-s(i  ,j,k  ))+&
                               dz(k-1,RIGHT_EDGE)*(s(i  ,j,k  )-s(i  ,j,k-1)))

              rPhiZW = (1./2)*(dz(k,  RIGHT_EDGE)*(s(i  ,j,k+1)-s(i  ,j,k  ))+&
                               dz(k-1,RIGHT_EDGE)*(s(i  ,j,k  )-s(i  ,j,k-1))+&
                               dz(k,  RIGHT_EDGE)*(s(i-1,j,k+1)-s(i-1,j,k  ))+&
                               dz(k-1,RIGHT_EDGE)*(s(i-1,j,k  )-s(i-1,j,k-1)))

              rPhiZN = (1./2)*(dz(k,  RIGHT_EDGE)*(s(i,j+1,k+1)-s(i,j+1,k  ))+&
                               dz(k-1,RIGHT_EDGE)*(s(i,j+1,k  )-s(i,j+1,k-1))+&
                               dz(k,  RIGHT_EDGE)*(s(i,j  ,k+1)-s(i,j  ,k  ))+&
                               dz(k-1,RIGHT_EDGE)*(s(i,j  ,k  )-s(i,j  ,k-1)))

              rPhiZS = (1./2)*(dz(k,  RIGHT_EDGE)*(s(i,j  ,k+1)-s(i,j  ,k  ))+&
                               dz(k-1,RIGHT_EDGE)*(s(i,j  ,k  )-s(i,j  ,k-1))+&
                               dz(k,  RIGHT_EDGE)*(s(i,j-1,k+1)-s(i,j-1,k  ))+&
                               dz(k-1,RIGHT_EDGE)*(s(i,j-1,k  )-s(i,j-1,k-1)))

              rPhiZF = dz(k,RIGHT_EDGE)*(s(i,j,k+1) - s(i,j,k)  )
              rPhiZB = dz(k,LEFT_EDGE)*(s(i,j,k)   - s(i,j,k-1))


              !- kpd - Compute the magnitude of the normal for ALL faces
              rMagE = sqrt( rPhiXE**2. + rPhiYE**2. + rPhiZE**2.) + eps
              rMagW = sqrt( rPhiXW**2. + rPhiYW**2. + rPhiZW**2.) + eps
              rMagN = sqrt( rPhiXN**2. + rPhiYN**2. + rPhiZN**2.) + eps
              rMagS = sqrt( rPhiXS**2. + rPhiYS**2. + rPhiZS**2.) + eps
              rMagF = sqrt( rPhiXF**2. + rPhiYF**2. + rPhiZF**2.) + eps
              rMagB = sqrt( rPhiXB**2. + rPhiYB**2. + rPhiZB**2.) + eps

              !------------------------------------------------------------
              !- kpd - Finally, compue the curvature, K=grad(s)/||grad(s)||
              crv(i,j,k) = dx(i,CENTER)*(rPhiXE/rMagE - rPhiXW/rMagW) &
                         + dy(j,CENTER)*(rPhiYN/rMagN - rPhiYS/rMagS) &
                         + dz(k,CENTER)*(rPhiZF/rMagF - rPhiZB/rMagB) 
              !------------------------------------------------------------

              end do
           end do
        end do

        !********************************************************************

        !----------------------------------------------------
        !- kpd - Set phase function on each side of interface
        !----------------------------------------------------
        do k = kz1-1,kz2+1
           do j = jy1-1,jy2+1
              do i = ix1-1,ix2+1
                 pf(i,j,k) = 0.

                 if(s(i,j,k).ge.0.) then
                    pf(i,j,k) = 1.                       
                    visc(i,j,k) = vis1/vis2               !- kpd - Set viscosity on each side of interface
                 else
                    visc(i,j,k) = vis2/vis2
                 end if

              end do
           end do
        end do

        !********************************************************************

        !--------------------------------------------------------------
        !- kpd - These are FACE VALUED inverse densities for each phase
        !--------------------------------------------------------------

        !- kpd - density on x-faces
        !rho1x = 0.
        !rho2x = 0.
        !- kpd - Loop through boundary and interior cell faces
        do k = kz1-1,kz2+1
           do j = jy1-1,jy2+1
              do i = ix1-1,ix2+1

              rho1x(i,j,k) = 0.
              rho2x(i,j,k) = 0.

              a1 = (pf(i-1,j,k) + pf(i,j,k)) / 2.                       
              a2 = pf(i-1,j,k)  /abs(pf(i-1,j,k)  +eps) * &
                   pf(i,j,k)/abs(pf(i,j,k)+eps)

              rho1x(i,j,k) = a1*a2/(rho1/rho2)
              rho2x(i,j,k) = (1. - a1*a2)/(rho2/rho2)

              end do
           end do
        end do

        !- kpd - density on y-faces
        !- kpd - Loop through boundary and interior cell faces
        do i = ix1-1,ix2+1
           do k = kz1-1,kz2+1
              do j = jy1-1,jy2+1

              rho1y(i,j,k) = 0.
              rho2y(i,j,k) = 0.

              a1 = (pf(i,j-1,k) + pf(i,j,k)) / 2.           
              a2 = pf(i,j-1,k)  /abs(pf(i,j-1,k)  +eps) * &
                   pf(i,j,k)/abs(pf(i,j,k)+eps)

              rho1y(i,j,k) = a1*a2/(rho1/rho2)
              rho2y(i,j,k) = (1. - a1*a2)/(rho2/rho2)

              end do
           end do
        end do

        !- kpd - density on z-faces
        !- kpd - Loop through boundary and interior cell faces
        do i = ix1-1,ix2+1
           do j = jy1-1,jy2+1
              do k = kz1-1,kz2+1

              rho1z(i,j,k) = 0.
              rho2z(i,j,k) = 0.

              a1 = (pf(i,j,k-1) + pf(i,j,k)) / 2.           
              a2 = pf(i,j,k-1)  /abs(pf(i,j,k-1)  +eps) * &
                   pf(i,j,k)/abs(pf(i,j,k)+eps)

              rho1z(i,j,k) = a1*a2/(rho1/rho2)
              rho2z(i,j,k) = (1. - a1*a2)/(rho2/rho2)

              end do
           end do
        end do

      end subroutine mph_KPDcurvature3DAB

!=========================================================================
!=========================================================================
!=========================================================================


!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------
!- kpd - There is a break in the ...curvature routine to 
!        go back and do guard cell filling.
!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------

!=========================================================================
!=========================================================================
!=========================================================================


        subroutine mph_KPDcurvature3DC(s,crv,rho1x,rho2x,rho1y,rho2y, &
                                       pf,w,wcd,sigx,sigy,sigcx,sigcy,dx,&
                                       dy,rho1,rho2,xit,ix1,ix2, &
                                       jy1,jy2,dz,kz1,kz2,rho1z, &
                                       rho2z,sigz,sigcz)


        implicit none

#include "Flash.h"
#include "constants.h"

        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        real, intent(in) :: rho1, rho2, xit

        real, dimension(:,:), intent(in) :: dx,dy,dz

        real, dimension(:,:,:), intent(inout):: s,crv,rho1x,rho2x,rho1y, &
                                                rho2y,pf,w,wcd,sigx,sigy, &
                                                sigcx,sigcy,rho1z,rho2z,sigz,sigcz

        !integer :: icrv(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
        integer :: icrv(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD)

        !- kpd - Local Data 
        real :: th,aa,xijl,xijr, &
                a1,a2,cri,xid,xij,xidl,xidr,yid,yidr,yidl,yij,yijl,yijr, &
                zijl,zijr,zid,zij,zidl,zidr
        integer :: i,j,k
        real :: dxci,dxci2,dycj,dycj2,dzck,dzck2,dxli,dylj,dzlk
        real, parameter :: eps = 1.0E-13


        sigx = 0.
        sigcx = 0.
        sigy = 0.
        sigcy = 0.
        sigz = 0.
        sigcz = 0.
        w = 0.
        wcd = 0.
        icrv = 0

        do k = kz1-1,kz2
           do j = jy1-1,jy2
              do i = ix1-1,ix2

              !--------------------------------------------------------------
              !- kpd - pf=0 in current cell and pf=1 in cell to right
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.0..and.pf(i+1,j,k).eq.1.) then

                 th = abs(s(i+1,j,k))/(abs(s(i+1,j,k))+abs(s(i,j,k)))

                 !- kpd - Unused in FLASH, needs to be phased out
                 cri = crv(i+1,j,k)*(1.-th) + crv(i,j,k)*th

                 xijl = xit*crv(i,j,k)                    !- kpd - sigma*K. Used for jump in pressure
                 xijr = xit*crv(i+1,j,k)                  !- kpd - sigma*K. Used for jump in pressure
                 xidl = 0.                                !- kpd - Used for jump in gradient
                 xidr = 0.                                !- kpd - Used for jump in gradient

                 xij = xijl*th + xijr*(1.-th)             !- kpd - Jump in value
                 xid = xidl*th + xidr*(1.-th)             !- kpd - Jump in gradient. Equal to 0 here.

                 !--------------------------------------------------------------------------
                 !- kpd - Density ALWAYS comes in as rho1x=0 and rho2x=1/rho2
                 !-----------------------------------------------------------
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)              !- kpd - Mixture density (Smeared)
                 rho1x(i+1,j,k) = rho1x(i+1,j,k)*(rho1/rho2)/aa  !- kpd - Density IS SMEARED HERE
                 rho2x(i+1,j,k) = rho2x(i+1,j,k)*(rho2/rho2)/aa  !- kpd - Density IS SMEARED HERE
                 !--------------------------------------------------------------------------

                 dxci = (1.0/dx(i,RIGHT_EDGE))*(1.0/dx(i,CENTER))
                 dxci2 =(1.0/dx(i+1,LEFT_EDGE))*(1.0/dx(i+1,CENTER))
                 dxli = 1.0/dx(i+1,LEFT_EDGE)

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k) - xij/aa/dxci !- xid*th*(rho1/rho2)/aa/dx
                 w(i+1,j,k) = w(i+1,j,k) + xij/aa/dxci2 !- xid*(1.-th)*(rho2/rho2)/aa/dx
                 wcd(i,j,k)   = wcd(i,j,k) - xij/dxci !- xid*th*(rho1/rho2)/aa/dx
                 wcd(i+1,j,k) = wcd(i+1,j,k) + xij/dxci2 !- xid*(1.-th)*(rho2/rho2)/aa/dx

                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigx(i+1,j,k) = -xij/aa/dxli
                 sigcx(i+1,j,k) = -xij/dxli


                 icrv(i,j,k) = 1
                 icrv(i+1,j,k) = 1

              end if

              !--------------------------------------------------------------
              !- kpd - pf=1 in current cell and pf=0 in cell to right
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.1..and.pf(i+1,j,k).eq.0.) then

                 th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k)))

                 !- kpd - Unused in FLASH, needs to be phased out
                 cri = crv(i,j,k)*(1.-th) + crv(i+1,j,k)*th

                 xijl = xit*crv(i,j,k)
                 xijr = xit*crv(i+1,j,k)
                 xidl = 0.
                 xidr = 0.

                 xij = xijl*(1.-th) + xijr*th
                 xid = xidl*(1.-th) + xidr*th

                 !--------------------------------------------------------------------------
                 !- kpd - Density ALWAYS comes in as rho1x=0 and rho2x=1/rho2
                 !-----------------------------------------------------------
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)              !- kpd - Density IS SMEARED HERE
                 rho1x(i+1,j,k) = rho1x(i+1,j,k)*(rho1/rho2)/aa  !- kpd - Density IS SMEARED HERE
                 rho2x(i+1,j,k) = rho2x(i+1,j,k)*(rho2/rho2)/aa  !- kpd - Density IS SMEARED HERE
                 !--------------------------------------------------------------------------
                 dxci = (1.0/dx(i,RIGHT_EDGE))*(1.0/dx(i,CENTER))
                 dxci2 =(1.0/dx(i+1,LEFT_EDGE))*(1.0/dx(i+1,CENTER))
                 dxli = 1.0/dx(i+1,LEFT_EDGE)

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k)   + xij/aa/dxci !+ xid*(1.-th)*(rho2/rho2)/aa/(1/dx(i,CENTER))
                 w(i+1,j,k) = w(i+1,j,k) - xij/aa/dxci2 !+ xid*th*(rho1/rho2)/aa/(1/dx(i+1,CENTER))
                 wcd(i,j,k)   = wcd(i,j,k)   + xij/dxci !+ xid*(1.-th)*(rho2/rho2)/aa/(1/dx(i,CENTER))
                 wcd(i+1,j,k) = wcd(i+1,j,k) - xij/dxci2 !+ xid*th*(rho1/rho2)/aa/(1/dx(i+1,CENTER))

                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigx(i+1,j,k) = xij/aa/dxli
                 sigcx(i+1,j,k) = xij/dxli

                 icrv(i,j,k) = 1
                 icrv(i+1,j,k) = 1

              end if

              !--------------------------------------------------------------
              !- kpd - pf=0 in current cell and pf=1 in cell above
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.0..and.pf(i,j+1,k).eq.1.) then

                 th = abs(s(i,j+1,k))/(abs(s(i,j+1,k))+abs(s(i,j,k)))

                 !- kpd - Unused in FLASH, needs to be phased out
                 cri = crv(i,j+1,k)*(1.-th) + crv(i,j,k)*th

                 yijl = xit*crv(i,j,k)
                 yijr = xit*crv(i,j+1,k)
                 yidl = 0.
                 yidr = 0.

                 yij = yijl*th + yijr*(1.-th)
                 yid = yidl*th + yidr*(1.-th)

                 !--------------------------------------------------------------------------
                 !- kpd - Density ALWAYS comes in as rho1y=0 and rho2y=1/rho2
                 !-----------------------------------------------------------
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)               !- kpd - Density IS SMEARED HERE
                 rho1y(i,j+1,k) = rho1y(i,j+1,k)*(rho1/rho2)/aa   !- kpd - Density IS SMEARED HERE
                 rho2y(i,j+1,k) = rho2y(i,j+1,k)*(rho2/rho2)/aa   !- kpd - Density IS SMEARED HERE

                 dycj = (1.0/dy(j,RIGHT_EDGE))*(1.0/dy(j,CENTER))
                 dycj2 =(1.0/dy(j+1,LEFT_EDGE))*(1.0/dy(j+1,CENTER))
                 dylj = 1.0/dy(j+1,LEFT_EDGE)

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k) - yij/aa/dycj !- yid*th*(rho1/rho2)/aa/(1/dy(j,CENTER))
                 w(i,j+1,k) = w(i,j+1,k)   + yij/aa/dycj2 !- yid*(1.-th)*(rho2/rho2)/aa/(1/dy(j+1,CENTER))
                 wcd(i,j,k)   = wcd(i,j,k) - yij/dycj !- yid*th*(rho1/rho2)/aa/(1/dy(j,CENTER))
                 wcd(i,j+1,k) = wcd(i,j+1,k)   + yij/dycj2 !- yid*(1.-th)*(rho2/rho2)/aa/(1/dy(j+1,CENTER))

                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigy(i,j+1,k) = - yij/aa/dylj
                 sigcy(i,j+1,k) = - yij/dylj

                 icrv(i,j,k) = 1
                 icrv(i,j+1,k) = 1

              end if

              !--------------------------------------------------------------
              !- kpd - pf=1 in current cell and pf=0 in cell above
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.1..and.pf(i,j+1,k).eq.0.) then

                 th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k)))

                 !- kpd - Unused in FLASH, needs to be phased out
                 cri = crv(i,j,k)*(1.-th) + crv(i,j+1,k)*th

                 yijl = xit*crv(i,j,k)
                 yijr = xit*crv(i,j+1,k)
                 yidl = 0.
                 yidr = 0.

                 yij = yijl*(1.-th) + yijr*th
                 yid = yidl*(1.-th) + yidr*th

                 !--------------------------------------------------------------------------
                 !- kpd - Density ALWAYS comes in as rho1y=0 and rho2y=1/rho2
                 !-----------------------------------------------------------
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)              !- kpd - Density IS SMEARED HERE
                 rho1y(i,j+1,k) = rho1y(i,j+1,k)*(rho1/rho2)/aa  !- kpd - Density IS SMEARED HERE
                 rho2y(i,j+1,k) = rho2y(i,j+1,k)*(rho2/rho2)/aa  !- kpd - Density IS SMEARED HERE

                 dycj = (1.0/dy(j,RIGHT_EDGE))*(1.0/dy(j,CENTER))
                 dycj2 =(1.0/dy(j+1,LEFT_EDGE))*(1.0/dy(j+1,CENTER))
                 dylj = 1.0/dy(j+1,LEFT_EDGE)

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k)   + yij/aa/dycj !+ yid*(1.-th)*(rho2/rho2)/aa/(1/dy(j,CENTER))
                 w(i,j+1,k) = w(i,j+1,k) - yij/aa/dycj2 !+ yid*th*(rho1/rho2)/aa/(1/dy(j+1,CENTER))
                 wcd(i,j,k)   = wcd(i,j,k)   + yij/dycj !+ yid*(1.-th)*(rho2/rho2)/aa/(1/dy(j,CENTER))
                 wcd(i,j+1,k) = wcd(i,j+1,k) - yij/dycj2 !+ yid*th*(rho1/rho2)/aa/(1/dy(j+1,CENTER))

                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigy(i,j+1,k) = yij/aa/dylj
                 sigcy(i,j+1,k) = yij/dylj

                 icrv(i,j,k) = 1
                 icrv(i,j+1,k) = 1

              end if

              !--------------------------------------------------------------
              !- kpd - pf=0 in current cell and pf=1 in cell to front
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.0..and.pf(i,j,k+1).eq.1.) then

                 th = abs(s(i,j,k+1))/(abs(s(i,j,k+1))+abs(s(i,j,k)))

                 !- kpd - Unused in FLASH, needs to be phased out
                 cri = crv(i,j,k+1)*(1.-th) + crv(i,j,k)*th

                 zijl = xit*crv(i,j,k)                    !- kpd - sigma*K. Used for jump in pressure
                 zijr = xit*crv(i,j,k+1)                  !- kpd - sigma*K. Used for jump in pressure
                 zidl = 0.                                !- kpd - Used for jump in gradient
                 zidr = 0.                                !- kpd - Used for jump in gradient

                 zij = zijl*th + zijr*(1.-th)             !- kpd - Jump in value
                 zid = zidl*th + zidr*(1.-th)             !- kpd - Jump in gradient. Equal to 0 here.

                 !--------------------------------------------------------------------------
                 !- kpd - Density ALWAYS comes in as rho1z=0 and rho2z=1/rho2
                 !-----------------------------------------------------------
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)              !- kpd - Mixture density (Smeared)
                 rho1z(i,j,k+1) = rho1z(i,j,k+1)*(rho1/rho2)/aa  !- kpd - Density IS SMEARED HERE  
                 rho2z(i,j,k+1) = rho2z(i,j,k+1)*(rho2/rho2)/aa  !- kpd - Density IS SMEARED HERE

                 dzck = (1.0/dz(k,RIGHT_EDGE))*(1.0/dz(k,CENTER))
                 dzck2 =(1.0/dz(k+1,LEFT_EDGE))*(1.0/dz(k+1,CENTER))
                 dzlk = 1.0/dz(k+1,LEFT_EDGE)

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k)   - zij/aa/dzck !- zid*th*(rho1/rho2)/aa/(1/dz(k,CENTER))
                 w(i,j,k+1) = w(i,j,k+1) + zij/aa/dzck2 !- zid*(1.-th)*(rho2/rho2)/aa/(1/dz(k+1,CENTER))
                 wcd(i,j,k)   = wcd(i,j,k)   - zij/dzck !- zid*th*(rho1/rho2)/aa/(1/dz(k,CENTER))
                 wcd(i,j,k+1) = wcd(i,j,k+1) + zij/dzck2 !- zid*(1.-th)*(rho2/rho2)/aa/(1/dz(k+1,CENTER))

                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigz(i,j,k+1) = - zij/aa/dzlk
                 sigcz(i,j,k+1) = - zij/dzlk

                 icrv(i,j,k)   = 1
                 icrv(i,j,k+1) = 1

              end if

              !--------------------------------------------------------------
              !- kpd - pf=1 in current cell and pf=0 in cell to front
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.1..and.pf(i,j,k+1).eq.0.) then

                 th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j,k+1)))

                 !- kpd - Unused in FLASH, needs to be phased out
                 cri = crv(i,j,k)*(1.-th) + crv(i,j,k+1)*th

                 zijl = xit*crv(i,j,k)
                 zijr = xit*crv(i,j,k+1)
                 zidl = 0.
                 zidr = 0.

                 zij = zijl*(1.-th) + zijr*th
                 zid = zidl*(1.-th) + zidr*th

                 !--------------------------------------------------------------------------
                 !- kpd - Density ALWAYS comes in as rho1z=0 and rho2z=1/rho2
                 !-----------------------------------------------------------
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)              !- kpd - Density IS SMEARED HERE
                 rho1z(i,j,k+1) = rho1z(i,j,k+1)*(rho1/rho2)/aa  !- kpd - Density IS SMEARED HERE
                 rho2z(i,j,k+1) = rho2z(i,j,k+1)*(rho2/rho2)/aa  !- kpd - Density IS SMEARED HERE

                 dzck = (1.0/dz(k,RIGHT_EDGE))*(1.0/dz(k,CENTER))
                 dzck2 =(1.0/dz(k+1,LEFT_EDGE))*(1.0/dz(k+1,CENTER))
                 dzlk = 1.0/dz(k+1,LEFT_EDGE)

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k)   + zij/aa/dzck !+ zid*(1.-th)*(rho2/rho2)/aa/(1/dz(k,CENTER))
                 w(i,j,k+1) = w(i,j,k+1) - zij/aa/dzck2 !+ zid*th*(rho1/rho2)/aa/(1/dz(k+1,CENTER))
                 wcd(i,j,k)   = wcd(i,j,k)   + zij/dzck !+ zid*(1.-th)*(rho2/rho2)/aa/(1/dz(k,CENTER))
                 wcd(i,j,k+1) = wcd(i,j,k+1) - zij/dzck2 !+ zid*th*(rho1/rho2)/aa/(1/dz(k+1,CENTER))

                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigz(i,j,k+1) = zij/aa/dzlk
                 sigcz(i,j,k+1) = zij/dzlk

                 icrv(i,j,k) = 1
                 icrv(i,j,k+1) = 1

              end if

              !--------------------------------------------------------------
              !--------------------------------------------------------------

           end do
        end do
        end do

        crv = crv*icrv

      end subroutine mph_KPDcurvature3DC

!=========================================================================
!=========================================================================
!=========================================================================

