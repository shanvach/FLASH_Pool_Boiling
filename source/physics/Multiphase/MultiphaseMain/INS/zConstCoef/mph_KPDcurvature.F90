!=========================================================================
!=========================================================================
!=========================================================================


        subroutine mph_KPDcurvature2DAB(s,crv,rho1x,rho2x,rho1y,rho2y,pf,w,sigx,sigy,dx,dy, &
           rho1,rho2,xit,crmx,crmn,ix1,ix2,jy1,jy2,visc,vis1,vis2)

   
        implicit none

#include "Flash.h"
#include "constants.h"

        !*************************************************************************

        !---------------------------------
        !- kpd - Data from routine call...
        !---------------------------------
        integer, intent(in) :: ix1,ix2,jy1,jy2
        real, intent(in) :: rho1, rho2, xit, vis1, vis2
        real, dimension(:,:), intent(in) :: dx,dy 
        real, intent(out) :: crmx, crmn
        real, dimension(:,:,:), intent(inout):: s,crv,rho1x,rho2x,rho1y, &
                                               rho2y,pf,w,sigx,sigy,visc

        !--------------------------
        !- kpd - Local variables...
        !--------------------------
        integer :: i,j,k
        real :: rPhiXN,rPhiXE,rPhiXS,rPhiXW, &
                rPhiYN,rPhiYE,rPhiYS,rPhiYW, &
                rMagN,rMagE,rMagS,rMagW

        real :: th,aa,xijl,xijr, &
                a1,a2,cri,xid,xij,xidl,xidr,yid,yidr,yidl,yij,yijl,yijr
        real, parameter :: eps = 1E-13

        real, dimension(ix2,jy2) :: adf
        real :: axr,axl,ayr,ayl,ax,ay

        real :: dsdx_ll,dsdx_lr,dsdx_rl,dsdx_rr,&
                rPhiXN_1,rPhiXN_2,rPhiXS_1,rPhiXS_2
        real :: dsdy_ll,dsdy_lr,dsdy_rl,dsdy_rr,&
                rPhiYE_1,rPhiYE_2,rPhiYW_1,rPhiYW_2
        real :: dxr,dxl,dyr,dyl

        crmx = -1E10
        crmn = 1E10


        !*************************************************************************

        !- kpd - For 2D runs 
        k=1

        !- kpd - Compute the curvature ---------------

        crv = 0.
        do j = jy1,jy2
           do i = ix1,ix2    
              !----------------------------------------------------
              !- kpd - 2 phi gradients per face method
              !----------------------------------------------------
              !        X - Location
              rPhiXE = dx(i,RIGHT_EDGE)*(s(i+1,j,k)-s(i,j,k)  )
              rPhiXW = dx(i,LEFT_EDGE)*(s(i,j,k)  -s(i-1,j,k))

              dsdx_rr = dx(i  ,RIGHT_EDGE)*(s(i+1,j+1,k)-s(i  ,j+1,k))
              dsdx_lr = dx(i  ,LEFT_EDGE )*(s(i  ,j+1,k)-s(i-1,j+1,k))
              dsdx_rl = dx(i  ,RIGHT_EDGE)*(s(i+1,j  ,k)-s(i  ,j  ,k))
              dsdx_ll = dx(i  ,LEFT_EDGE )*(s(i  ,j  ,k)-s(i-1,j  ,k))

              dxr = (1/dx(i  ,CENTER))/2
              dyr = (1/dy(j+1,CENTER))/2
              dyl = (1/dy(j  ,CENTER))/2

              rPhiXN_1 = (dxr/(2*dxr))*dsdx_ll + (dxr/(2*dxr))*dsdx_rl
              rPhiXN_2 = (dxr/(2*dxr))*dsdx_lr + (dxr/(2*dxr))*dsdx_rr

              ! Bilinear Interpolation
              rPhiXN   = (dyr/(dyl+dyr))*rPhiXN_1 + (dyl/(dyl+dyr))*rPhiXN_2

              ! Simple Average
              !rPhiXN = (1./4)*(dx(i  ,RIGHT_EDGE)*(s(i+1,j+1,k)-s(i  ,j+1,k))+&
              !                 dx(i  ,LEFT_EDGE )*(s(i  ,j+1,k)-s(i-1,j+1,k))+&
              !                 dx(i  ,RIGHT_EDGE)*(s(i+1,j  ,k)-s(i  ,j  ,k))+&
              !                 dx(i  ,LEFT_EDGE )*(s(i  ,j  ,k)-s(i-1,j  ,k)))

              dsdx_rr = dx(i,  RIGHT_EDGE)*(s(i+1,j  ,k)-s(i  ,j  ,k))
              dsdx_lr = dx(i  ,LEFT_EDGE )*(s(i  ,j  ,k)-s(i-1,j  ,k))
              dsdx_rl = dx(i,  RIGHT_EDGE)*(s(i+1,j-1,k)-s(i  ,j-1,k))
              dsdx_ll = dx(i  ,LEFT_EDGE )*(s(i  ,j-1,k)-s(i-1,j-1,k))

              dxr = (1/dx(i  ,CENTER))/2
              dyr = (1/dy(j  ,CENTER))/2
              dyl = (1/dy(j-1,CENTER))/2

              rPhiXS_1 = (dxr/(2*dxr))*dsdx_ll + (dxr/(2*dxr))*dsdx_rl
              rPhiXS_2 = (dxr/(2*dxr))*dsdx_lr + (dxr/(2*dxr))*dsdx_rr

              ! Bilinear Interpolation
              rPhiXS   = (dyr/(dyl+dyr))*rPhiXS_1 + (dyl/(dyl+dyr))*rPhiXS_2

              ! Simple Average
              !rPhiXS = (1./4)*(dx(i,  RIGHT_EDGE)*(s(i+1,j  ,k)-s(i  ,j  ,k))+&
              !                 dx(i  ,LEFT_EDGE )*(s(i  ,j  ,k)-s(i-1,j  ,k))+&
              !                 dx(i,  RIGHT_EDGE)*(s(i+1,j-1,k)-s(i  ,j-1,k))+&
              !                 dx(i  ,LEFT_EDGE )*(s(i  ,j-1,k)-s(i-1,j-1,k)))

             !        Y - Location
              rPhiYN = dy(j,RIGHT_EDGE)*(s(i,j+1,k)-s(i,j,k)  )
              rPhiYS = dy(j,LEFT_EDGE)*(s(i,j,k)  -s(i,j-1,k))

              dsdy_rr = dy(j  ,RIGHT_EDGE)*(s(i+1,j+1,k)-s(i+1,j  ,k))
              dsdy_rl = dy(j  ,LEFT_EDGE )*(s(i+1,j  ,k)-s(i+1,j-1,k))
              dsdy_lr = dy(j  ,RIGHT_EDGE)*(s(i  ,j+1,k)-s(i  ,j  ,k))
              dsdy_ll = dy(j  ,LEFT_EDGE )*(s(i  ,j  ,k)-s(i  ,j-1,k))

              dxr = (1/dx(i+1,CENTER))/2
              dxl = (1/dx(i  ,CENTER))/2
              dyr = (1/dy(j  ,CENTER))/2

              rPhiYE_1 = (dyr/(2*dyr))*dsdy_ll + (dyr/(2*dyr))*dsdy_lr
              rPhiYE_2 = (dyr/(2*dyr))*dsdy_rl + (dyr/(2*dyr))*dsdy_rr

              ! Bilinear Interpolation
              rPhiYE   = (dxr/(dxr+dxl))*rPhiYE_1 + (dxl/(dxr+dxl))*rPhiYE_2

              ! Simple Average
              !rPhiYE = (1./4)*(dy(j  ,RIGHT_EDGE)*(s(i+1,j+1,k)-s(i+1,j  ,k))+&
              !                 dy(j  ,LEFT_EDGE )*(s(i+1,j  ,k)-s(i+1,j-1,k))+&
              !                 dy(j  ,RIGHT_EDGE)*(s(i  ,j+1,k)-s(i  ,j  ,k))+&
              !                 dy(j  ,LEFT_EDGE )*(s(i  ,j  ,k)-s(i  ,j-1,k)))

              dsdy_rr = dy(j  ,RIGHT_EDGE)*(s(i  ,j+1,k)-s(i  ,j  ,k))
              dsdy_rl = dy(j  ,LEFT_EDGE )*(s(i  ,j  ,k)-s(i  ,j-1,k))
              dsdy_lr = dy(j  ,RIGHT_EDGE)*(s(i-1,j+1,k)-s(i-1,j  ,k))
              dsdy_ll = dy(j  ,LEFT_EDGE )*(s(i-1,j  ,k)-s(i-1,j-1,k))

              dxr = (1/dx(i  ,CENTER))/2
              dxl = (1/dx(i-1,CENTER))/2
              dyr = (1/dy(j  ,CENTER))/2

              rPhiYW_1 = (dyr/(2*dyr))*dsdy_ll + (dyr/(2*dyr))*dsdy_lr
              rPhiYW_2 = (dyr/(2*dyr))*dsdy_rl + (dyr/(2*dyr))*dsdy_rr

              ! Bilinear Interpolation
              rPhiYW   = (dxr/(dxr+dxl))*rPhiYW_1 + (dxl/(dxr+dxl))*rPhiYW_2

              ! Simple Average
              !rPhiYW = (1./4)*(dy(j  ,RIGHT_EDGE)*(s(i  ,j+1,k)-s(i  ,j  ,k))+&
              !                 dy(j  ,LEFT_EDGE )*(s(i  ,j  ,k)-s(i  ,j-1,k))+&
              !                 dy(j  ,RIGHT_EDGE)*(s(i-1,j+1,k)-s(i-1,j  ,k))+&
              !                 dy(j  ,LEFT_EDGE )*(s(i-1,j  ,k)-s(i-1,j-1,k)))

             !----------------------------------------------------

              !- kpd - Compute the magnitude of the gradient at each face
              rMagE = sqrt( rPhiXE**2. + rPhiYE**2. ) + eps
              rMagW = sqrt( rPhiXW**2. + rPhiYW**2. ) + eps
              rMagN = sqrt( rPhiXN**2. + rPhiYN**2. ) + eps
              rMagS = sqrt( rPhiXS**2. + rPhiYS**2. ) + eps


              crv(i,j,k) = dx(i,CENTER) * (rPhiXE/rMagE - rPhiXW/rMagW) &
                         + dy(j,CENTER) * (rPhiYN/rMagN - rPhiYS/rMagS)
              !----------------------------------------------------

           end do
        end do

!       !---- RIAZ'S IMPLEMENTATION ----!
!
!        do j = jy1,jy2
!           do i = ix1,ix2
!
!              adf(i,j) = sqrt( ((s(i+1,j,k)-s(i-1,j,k))/2./dx)**2 + &
!                               ((s(i,j+1,k)-s(i,j-1,k))/2./dy)**2 )
!
!           end do
!        end do
!
!        crv = 0.
!        do j = jy1+1,jy2-1
!           do i = ix1+1,ix2-1
!
!              axr = (adf(i,j)+adf(i+1,j))/2.
!              axl = (adf(i-1,j)+adf(i,j))/2.
!              ayr = (adf(i,j)+adf(i,j+1))/2.
!              ayl = (adf(i,j-1)+adf(i,j))/2.
!
!              ax = (s(i+1,j,k)-s(i,j,k))/axr - (s(i,j,k)-s(i-1,j,k))/axl
!              ay = (s(i,j+1,k)-s(i,j,k))/ayr - (s(i,j,k)-s(i,j-1,k))/ayl
!
!              crv(i,j,k) = ax/dx**2 + ay/dy**2
!
!           end do
!        end do

        !*************************************************************************

        !----------------------------------
        !- kpd - Compute the phase function
        !----------------------------------
        k=1
        do j = jy1-1,jy2+1
           do i = ix1-1,ix2+1
              pf(i,j,k) = 0.

              if(s(i,j,k).ge.0.) then
                 pf(i,j,k) = 1.                       !- kpd - Set phase function on each side of interface
                 visc(i,j,k) = vis1/vis2               !- kpd - Set viscosity on each side of interface
              else
                 visc(i,j,k) = vis2/vis2
              end if
           end do
        end do

        !--------------------------------------------------------------
        !- kpd - These are FACE VALUED inverse densities for each phase
        !--------------------------------------------------------------

        !- kpd - density on x-face
        rho1x = 0.
        rho2x = 0.
        !- kpd - Loop through boundary and interior cell faces
        do j = jy1-1,jy2+1
           do i = ix1-1,ix2+1
              a1 = (pf(i-1,j,k) + pf(i,j,k)) / 2.                       
              a2 = pf(i-1,j,k)  /abs(pf(i-1,j,k)  +eps) * &
                   pf(i,j,k)/abs(pf(i,j,k)+eps)
              rho1x(i,j,k) = a1*a2/(rho1/rho2)
              rho2x(i,j,k) = (1. - a1*a2)/(rho2/rho2)
           end do
        end do

        !- kpd - density on y-face
        rho1y = 0.
        rho2y = 0.
        !- kpd - Loop through boundary and interior cell faces
        do i = ix1-1,ix2+1
           do j = jy1-1,jy2+1
              a1 = (pf(i,j-1,k) + pf(i,j,k)) / 2.           
              a2 = pf(i,j-1,k)  /abs(pf(i,j-1,k)  +eps) * &
                   pf(i,j,k)/abs(pf(i,j,k)+eps)
              rho1y(i,j,k) = a1*a2/(rho1/rho2)
              rho2y(i,j,k) = (1. - a1*a2)/(rho2/rho2)
           end do
        end do

      end subroutine mph_KPDcurvature2DAB

!=========================================================================
!=========================================================================
!=========================================================================
!=========================================================================
!=========================================================================
!=========================================================================

!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------
!- kpd - There is a break in the ...curvature routine to 
!        go back and do guard cell filling for curvature.
!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------

!=========================================================================
!=========================================================================
!=========================================================================
!- kpd - From here on mph_curvature2DC:
!=========================================================================
!=========================================================================
!=========================================================================

        subroutine mph_KPDcurvature2DC(s,crv,rho1x,rho2x,rho1y,rho2y,pf,w,wcd,sigx,sigy,&
                                sigcx,sigcy,dx,dy,rho1,rho2,xit,crmx,crmn,ix1,ix2,jy1,jy2)   

   
        use Multiphase_data, ONLY : mph_meshMe

        implicit none

#include "Flash.h"
#include "constants.h"

        !implicit real*8(a-h,o-z)
        !common/param/re,xit,rho1,rho2,g,sp,ubc,uout,nint

        integer, intent(in) :: ix1,ix2,jy1,jy2
        real, intent(in) :: rho1, rho2, xit
        real, dimension(:,:), intent(in) :: dx,dy 
        real, intent(out) :: crmx, crmn

        real, dimension(:,:,:), intent(inout):: s,crv,rho1x,rho2x,rho1y, &
                                               rho2y,pf,w,wcd,sigx,sigy,sigcx,sigcy
        integer :: icrv(NXB+2*NGUARD,NYB+2*NGUARD,1)

        !- kpd - 
        real :: th,aa,xijl,xijr, &
                a1,a2,cri,xid,xij,xidl,xidr,yid,yidr,yidl,yij,yijl,yijr
        integer :: i,j,k
        real, parameter :: eps = 1E-13

        integer :: iSmear

        real :: dxci,dxci2,dycj,dycj2,dxli,dylj

!--------------------------------------------
!----------------jump conditions ------------
!--------------------------------------------
!xij: jump in value
!xid: jump in gradient 
!l,r, values at pts left and right of the interface 
!crv: curvature
!cri: curvature at interface 
!left interface between i and i+1
!  phase 2 at i, phase 1 at i+1
!  theta = x_i+1 - x_I
!right interface between i and i+1
!  phase 1 at i, phase 2 at i+1
!  theta = x_I - x_i  
!--------------------------------------------
!--------------------------------------------

        iSmear  = 1

        crmx = -1E10
        crmn = 1E10
        sigx = 0.
        sigy = 0.
        sigcx = 0.
        sigcy = 0.
        w = 0.
        wcd = 0.
        icrv = 0

        !- kpd - Need to loop through one guard cell on each side to set jumps 
        !           when they cross block boundaries
                    !do j = jy1,jy2
                    !   do i = ix1,ix2
        k=1
        do j = jy1-1,jy2
           do i = ix1-1,ix2

              !--------------------------------------------------------------
              !- kpd - pf=0 (water) in current cell and pf=1 (air) in cell to right
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.0..and.pf(i+1,j,k).eq.1.) then

                 !          = (+)            = (+)           = (-)
                 th = abs(s(i+1,j,k))/(abs(s(i+1,j,k))+abs(s(i,j,k)))
                 
                 cri = crv(i+1,j,k)*(1.-th) + crv(i,j,k)*th

                 xijl = xit*crv(i,j,k)                 !- kpd - sigma*K. Used for jump in pressure
                 xijr = xit*crv(i+1,j,k)               !- kpd - sigma*K. Used for jump in pressure
                 xidl = 0.                             !- kpd - Used for jump in gradient
                 xidr = 0.                             !- kpd - Used for jump in gradient

                 xij = xijl*th + xijr*(1.-th)          !- kpd - Jump in value
                 xid = xidl*th + xidr*(1.-th)          !- kpd - Jump in gradient. Equal to 0 here.
                 
                 if (iSmear  .eq. 1) then

                 !- kpd - All Densities are relative to rho2...
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)           !- kpd - Mixture density (not inverse)

                 rho1x(i+1,j,k) = rho1x(i+1,j,k)*(rho1/rho2)/aa
                 rho2x(i+1,j,k) = rho2x(i+1,j,k)*(rho2/rho2)/aa 

                 dxci = (1.0/dx(i,RIGHT_EDGE))*(1.0/dx(i,CENTER))
                 dxci2 =(1.0/dx(i+1,LEFT_EDGE))*(1.0/dx(i+1,CENTER))
                 dxli = 1.0/dx(i+1,LEFT_EDGE)

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k)   - xij/aa/dxci - xid*th*(rho1/rho2)/aa/(1/dx(i,CENTER))
                 w(i+1,j,k) = w(i+1,j,k) + xij/aa/dxci2 - xid*(1.-th)*(rho2/rho2)/aa/(1/dx(i+1,CENTER))
                 wcd(i,j,k)   = wcd(i,j,k)   - xij/dxci - xid*th*(rho1/rho2)/(1/dx(i,CENTER))  
                 wcd(i+1,j,k) = wcd(i+1,j,k) + xij/dxci2 - xid*(1.-th)*(rho2/rho2)/(1/dx(i+1,CENTER))
                  
                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigx(i+1,j,k) = - xij/aa/dxli           !- kpd - sigma*K/rho/dx 
                 sigcx(i+1,j,k) = -xij/dxli

                 else

                 end if

                 crmx = max(abs(cri),crmx)
                 crmn = min(abs(cri),crmn)

                 icrv(i,j,k) = 1
                 icrv(i+1,j,k) = 1

              end if

              !--------------------------------------------------------------
              !- kpd - pf=1 in current cell and pf=0 in cell to right
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.1..and.pf(i+1,j,k).eq.0.) then

                 th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k)))

                 cri = crv(i,j,k)*(1.-th) + crv(i+1,j,k)*th

                 xijl = xit*crv(i,j,k)
                 xijr = xit*crv(i+1,j,k)
                 xidl = 0. 
                 xidr = 0.

                 xij = xijl*(1.-th) + xijr*th
                 xid = xidl*(1.-th) + xidr*th

                 if (iSmear  .eq. 1) then

                 !- kpd - All Densities are relative to rho2...
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)

                 rho1x(i+1,j,k) = rho1x(i+1,j,k)*(rho1/rho2)/aa
                 rho2x(i+1,j,k) = rho2x(i+1,j,k)*(rho2/rho2)/aa   

                 dxci = (1.0/dx(i,RIGHT_EDGE))*(1.0/dx(i,CENTER))
                 dxci2 =(1.0/dx(i+1,LEFT_EDGE))*(1.0/dx(i+1,CENTER))
                 dxli = 1.0/dx(i+1,LEFT_EDGE)

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k)   + xij/aa/dxci + xid*(1.-th)*(rho2/rho2)/aa/(1/dx(i,CENTER))
                 w(i+1,j,k) = w(i+1,j,k) - xij/aa/dxci2 + xid*th*(rho1/rho2)/aa/(1/dx(i+1,CENTER))
                 wcd(i,j,k)   = wcd(i,j,k)   + xij/dxci + xid*(1.-th)*(rho2/rho2)/(1/dx(i,CENTER))
                 wcd(i+1,j,k) = wcd(i+1,j,k) - xij/dxci2 + xid*th*(rho1/rho2)/(1/dx(i+1,CENTER))
                  
                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigx(i+1,j,k) = xij/aa/dxli
                 sigcx(i+1,j,k) = xij/dxli

                 else

                 end if

                 crmx = max(abs(cri),crmx)
                 crmn = min(abs(cri),crmn)

                 icrv(i,j,k) = 1
                 icrv(i+1,j,k) = 1

              end if

              !--------------------------------------------------------------
              !- kpd - pf=0 in current cell and pf=1 in cell above
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.0..and.pf(i,j+1,k).eq.1.) then

                 th = abs(s(i,j+1,k))/(abs(s(i,j+1,k))+abs(s(i,j,k)))

                 cri = crv(i,j+1,k)*(1.-th) + crv(i,j,k)*th

                 yijl = xit*crv(i,j,k)
                 yijr = xit*crv(i,j+1,k)
                 yidl = 0. 
                 yidr = 0.

                 yij = yijl*th + yijr*(1.-th)
                 yid = yidl*th + yidr*(1.-th)

                 if (iSmear  .eq. 1) then

                 !- kpd - All Densities are relative to rho2...
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)

                 rho1y(i,j+1,k) = rho1y(i,j+1,k)*(rho1/rho2)/aa
                 rho2y(i,j+1,k) = rho2y(i,j+1,k)*(rho2/rho2)/aa 

                 dycj = (1.0/dy(j,RIGHT_EDGE))*(1.0/dy(j,CENTER))
                 dycj2 =(1.0/dy(j+1,LEFT_EDGE))*(1.0/dy(j+1,CENTER))
                 dylj = 1.0/dy(j+1,LEFT_EDGE)

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k) - yij/aa/dycj - yid*th*(rho1/rho2)/aa/(1/dy(j,CENTER))
                 w(i,j+1,k) = w(i,j+1,k)   + yij/aa/dycj2 - yid*(1.-th)*(rho2/rho2)/aa/(1/dy(j+1,CENTER)) 
                 wcd(i,j,k)   = wcd(i,j,k) - yij/dycj - yid*th*(rho1/rho2)/(1/dy(j,CENTER))
                 wcd(i,j+1,k) = wcd(i,j+1,k) + yij/dycj2 - yid*(1.-th)*(rho2/rho2)/(1/dy(j+1,CENTER)) 
                  
                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigy(i,j+1,k) = - yij/aa/dylj
                 sigcy(i,j+1,k) = - yij/dylj

                 else

                 end if

                 crmx = max(abs(cri),crmx)
                 crmn = min(abs(cri),crmn)

                 icrv(i,j,k) = 1
                 icrv(i,j+1,k) = 1

              end if

              !--------------------------------------------------------------
              !- kpd - pf=1 in current cell and pf=0 in cell above
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.1..and.pf(i,j+1,k).eq.0.) then

                 th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k))) 
 
                 cri = crv(i,j,k)*(1.-th) + crv(i,j+1,k)*th

                 yijl = xit*crv(i,j,k)
                 yijr = xit*crv(i,j+1,k)
                 yidl = 0. 
                 yidr = 0.

                 yij = yijl*(1.-th) + yijr*th
                 yid = yidl*(1.-th) + yidr*th

                 if (iSmear  .eq. 1) then

                 !- kpd - All Densities are relative to rho2...
                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)

                 rho1y(i,j+1,k) = rho1y(i,j+1,k)*(rho1/rho2)/aa
                 rho2y(i,j+1,k) = rho2y(i,j+1,k)*(rho2/rho2)/aa  

                 dycj = (1.0/dy(j,RIGHT_EDGE))*(1.0/dy(j,CENTER))
                 dycj2 =(1.0/dy(j+1,LEFT_EDGE))*(1.0/dy(j+1,CENTER))
                 dylj = 1.0/dy(j+1,LEFT_EDGE)

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k)   + yij/aa/dycj + yid*(1.-th)*(rho2/rho2)/aa/(1/dy(j,CENTER))
                 w(i,j+1,k) = w(i,j+1,k) - yij/aa/dycj2 + yid*th*(rho1/rho2)/aa/(1/dy(j+1,CENTER))
                 wcd(i,j,k)   = wcd(i,j,k)   + yij/dycj + yid*(1.-th)*(rho2/rho2)/(1/dy(j,CENTER))
                 wcd(i,j+1,k) = wcd(i,j+1,k) - yij/dycj2 + yid*th*(rho1/rho2)/(1/dy(j+1,CENTER))
                  
                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigy(i,j+1,k) = yij/aa/dylj
                 sigcy(i,j+1,k) = yij/dylj

                 else

                 end if

                 crmx = max(abs(cri),crmx)
                 crmn = min(abs(cri),crmn)

                 icrv(i,j,k) = 1
                 icrv(i,j+1,k) = 1

              end if

              !--------------------------------------------------------------
              !--------------------------------------------------------------


           end do
        end do

        !- kpd - This is done for post-processing to visualize where jumps were applied...
        crv = crv*icrv

      end subroutine mph_KPDcurvature2DC


