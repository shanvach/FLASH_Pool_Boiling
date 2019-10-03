!=========================================================================
!=========================================================================
!=========================================================================


        subroutine mph_KPDcurvature2DAB(s,crv,rho1x,rho2x,rho1y,rho2y,pf,w,sigx,sigy,dx,dy, &
           rho1,rho2,xit,crmx,crmn,ix1,ix2,jy1,jy2,visc,vis1,vis2,alph,thco1,thco2,cp1,cp2,nrmx,nrmy,mflg,smhv,smrh,&
           lambda)

   
        implicit none

#include "Flash.h"
#include "constants.h"
        !*************************************************************************

        !---------------------------------
        !- kpd - Data from routine call...
        !---------------------------------
        integer, intent(in) :: ix1,ix2,jy1,jy2
        real, intent(in) :: dx,dy,rho1,rho2,xit,vis1,vis2,thco1,thco2,cp1,cp2 
        real, intent(out) :: crmx, crmn
        real, dimension(:,:,:), intent(inout):: s,crv,rho1x,rho2x,rho1y, &
                                               rho2y,pf,w,sigx,sigy,visc,&
                                               alph,nrmx,nrmy,mflg,smhv,smrh,lambda

        !--------------------------
        !- kpd - Local variables...
        !--------------------------
        integer :: i,j,k,kz1
        real :: rPhiXN,rPhiXE,rPhiXS,rPhiXW, &
                rPhiYN,rPhiYE,rPhiYS,rPhiYW, &
                rMagN,rMagE,rMagS,rMagW

        real :: th,aa,xijl,xijr, &
                a1,a2,cri,xid,xij,xidl,xidr,yid,yidr,yidl,yij,yijl,yijr
        real, parameter :: eps = 1E-13

        real, dimension(ix2,jy2) :: adf
        real :: axr,axl,ayr,ayl,ax,ay
        real :: sunion(NXB+2*NGUARD,NYB+2*NGUARD,1)
        real :: pfl(NXB+2*NGUARD,NYB+2*NGUARD,1)
        real :: b1,b2
        real :: rho3, vis3, thco3, cp3

        rho3  = (rho1 + rho2)/2
        vis3  = (vis1 + vis2)/2
        thco3 = thco2
        cp3   = cp2

        crmx = -1E10
        crmn = 1E10

        mflg = 0.0
        pfl  = 0.0

        !*************************************************************************

        !- kpd - For 2D runs 
        k=1
        kz1 = k

        sunion = s

        !do j=jy1-2,jy2+2
        !   do i=ix1-2,ix2+2
        !       sunion(i,j,kz1) = min(s(i,j,kz1),-lambda(i,j,kz1))
        !   end do
        !end do

        pf(ix1-1:ix2+1,jy1-1:jy2+1,k)   = 0.0
        pf(ix1-1:ix2+1,jy1-1:jy2+1,k)   = (sign(1.0,sunion(ix1-1:ix2+1,jy1-1:jy2+1,k))+1.0)/2.0

        pfl(ix1-1:ix2+1,jy1-1:jy2+1,k)  = 0.0
        pfl(ix1-1:ix2+1,jy1-1:jy2+1,k)  = (sign(1.0,lambda(ix1-1:ix2+1,jy1-1:jy2+1,k))+1.0)/2.0

        crv = 0.
        do j = jy1,jy2
           do i = ix1,ix2
              !----------------------------------------------------
              !- kpd - 2 phi gradients per face method
              !----------------------------------------------------
              !        X - Location
              rPhiXE = 1./dx*(s(i+1,j,k)-s(i,j,k)  )
              rPhiXW = 1./dx*(s(i,j,k)  -s(i-1,j,k))
              rPhiXN = 1./4./dx * ( (s(i+1,j+1,k) - s(i-1,j+1,k)) &
                                  + (s(i+1,j,k)   - s(i-1,j,k)  ) )
              rPhiXS = 1./4./dx * ( (s(i+1,j,k)   - s(i-1,j,k)  ) &
                                  + (s(i+1,j-1,k) - s(i-1,j-1,k)) )
             !        Y - Location
              rPhiYN = 1./dy*(s(i,j+1,k)-s(i,j,k)  )
              rPhiYS = 1./dy*(s(i,j,k)  -s(i,j-1,k))
              rPhiYE = 1./4./dy * ( (s(i+1,j+1,k) - s(i+1,j-1,k)) &
                                  + (s(i,j+1,k)   - s(i,j-1,k)  ) )
              rPhiYW = 1./4./dy * ( (s(i,j+1,k)   - s(i,j-1,k)  ) &
                                  + (s(i-1,j+1,k) - s(i-1,j-1,k)) )
              !----------------------------------------------------

              !- kpd - Compute the magnitude of the gradient at each face
              rMagE = sqrt( rPhiXE**2. + rPhiYE**2. ) + eps
              rMagW = sqrt( rPhiXW**2. + rPhiYW**2. ) + eps
              rMagN = sqrt( rPhiXN**2. + rPhiYN**2. ) + eps
              rMagS = sqrt( rPhiXS**2. + rPhiYS**2. ) + eps


              crv(i,j,k) = 1./dx * (rPhiXE/rMagE - rPhiXW/rMagW) &
                         + 1./dy * (rPhiYN/rMagN - rPhiYS/rMagS)
              !----------------------------------------------------

           end do
        end do

        do j=jy1-1,jy2+1
          do i=ix1-1,ix2+1

              visc(i,j,k) = (1-smhv(i,j,k))*(1-pfl(i,j,k))*(vis2/vis2) + &
                              (smhv(i,j,k))*(1-pfl(i,j,k))*(vis1/vis2) + &
                             (pfl(i,j,k))*(vis3/vis2)

              alph(i,j,k) = (1-pf(i,j,k))*(1-pfl(i,j,k))*(thco2/cp2)/(thco2/cp2) + &
                              (pf(i,j,k))*(1-pfl(i,j,k))*(thco1/cp1)/(thco2/cp2) + &
                             (pfl(i,j,k))*(thco3/cp3)/(thco2/cp2)

          end do
        end do

        do j = jy1+1,jy2-1
         do i = ix1+1,ix2-1

              if(s(i,j,k)*s(i+1,j,k) .le. 0.) mflg(i+1,j,k)   = 1.0
              if(s(i,j,k)*s(i,j+1,k) .le. 0.) mflg(i,j+1,k)   = 1.0
              if(s(i,j,k)*s(i-1,j,k) .le. 0.) mflg(i-1,j,k)   = 1.0
              if(s(i,j,k)*s(i,j-1,k) .le. 0.) mflg(i,j-1,k)   = 1.0

              if(s(i,j,k)*s(i+2,j,k) .le. 0.) mflg(i+2,j,k)   = 1.0
              if(s(i,j,k)*s(i,j+2,k) .le. 0.) mflg(i,j+2,k)   = 1.0
              if(s(i,j,k)*s(i-2,j,k) .le. 0.) mflg(i-2,j,k)   = 1.0
              if(s(i,j,k)*s(i,j-2,k) .le. 0.) mflg(i,j-2,k)   = 1.0

              if(s(i,j,k)*s(i+3,j,k) .le. 0.) mflg(i+3,j,k)   = 1.0
              if(s(i,j,k)*s(i,j+3,k) .le. 0.) mflg(i,j+3,k)   = 1.0
              if(s(i,j,k)*s(i-3,j,k) .le. 0.) mflg(i-3,j,k)   = 1.0
              if(s(i,j,k)*s(i,j-3,k) .le. 0.) mflg(i,j-3,k)   = 1.0

              if(s(i,j,k)*s(i+3,j,k) .le. 0.) mflg(i+4,j,k)   = 1.0
              if(s(i,j,k)*s(i,j+3,k) .le. 0.) mflg(i,j+4,k)   = 1.0
              if(s(i,j,k)*s(i-3,j,k) .le. 0.) mflg(i-4,j,k)   = 1.0
              if(s(i,j,k)*s(i,j-3,k) .le. 0.) mflg(i,j-4,k)   = 1.0

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

               b1 = 1. - (pfl(i-1,j,k) + pfl(i,j,k)) / 2.                      
               b2 = (1. - pfl(i-1,j,k)  /abs(pfl(i-1,j,k)  +eps)) * &
                    (1. - pfl(i,j,k)/abs(pfl(i,j,k)+eps))


               rho1x(i,j,k) = (a1*a2*(b1*b2))/(rho1/rho2)
               rho2x(i,j,k) = ((1. - a1*a2)*(b1*b2))/(rho2/rho2) + (1. - b1*b2)/(rho3/rho2)

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

              b1 = 1. - (pfl(i,j-1,k) + pfl(i,j,k)) / 2.           
              b2 = (1. - pfl(i,j-1,k)  /abs(pfl(i,j-1,k)  +eps)) * &
                   (1. - pfl(i,j,k)/abs(pfl(i,j,k)+eps))

              rho1y(i,j,k) = (a1*a2*(b1*b2))/(rho1/rho2)
              rho2y(i,j,k) = ((1. - a1*a2)*(b1*b2))/(rho2/rho2) + (1. - b1*b2)/(rho3/rho2)

           end do
        end do

        nrmx(ix1:ix2,jy1:jy2,kz1) =           &
           (( s(ix1+1:ix2+1,jy1:jy2,kz1) -   &
             s(ix1-1:ix2-1,jy1:jy2,kz1) )/2./dx)/ &
             sqrt( ((s(ix1+1:ix2+1,jy1:jy2,kz1) - &
             s(ix1-1:ix2-1,jy1:jy2,kz1))/2./dx)**2 &
             + ((s(ix1:ix2,jy1+1:jy2+1,kz1) - &
             s(ix1:ix2,jy1-1:jy2-1,kz1))/2./dy)**2 )
        
        nrmy(ix1:ix2,jy1:jy2,kz1) =           &
           (( s(ix1:ix2,jy1+1:jy2+1,kz1) -   &
             s(ix1:ix2,jy1-1:jy2-1,kz1) )/2./dy)/ &
             sqrt( ((s(ix1+1:ix2+1,jy1:jy2,kz1) - &
             s(ix1-1:ix2-1,jy1:jy2,kz1))/2./dx)**2 &
             + ((s(ix1:ix2,jy1+1:jy2+1,kz1) - &
             s(ix1:ix2,jy1-1:jy2-1,kz1))/2./dy)**2 )

        !nrmx(ix1-1:ix2+1,jy1-1:jy2+1,kz1) =           &
        !   (( s(ix1:ix2+2,jy1-1:jy2+1,kz1) -   &
        !      s(ix1-2:ix2,jy1-1:jy2+1,kz1) )/2./dx)/ &
        !     sqrt( ((s(ix1:ix2+2,jy1-1:jy2+1,kz1) - &
        !             s(ix1-2:ix2,jy1-1:jy2+1,kz1))/2./dx)**2 &
        !         + ((s(ix1-1:ix2+1,jy1:jy2+2,kz1) - &
        !             s(ix1-1:ix2+1,jy1-2:jy2,kz1))/2./dy)**2 )

        !nrmy(ix1-1:ix2+1,jy1-1:jy2+1,kz1) =           &
        !   (( s(ix1-1:ix2+1,jy1:jy2+2,kz1) -   &
        !      s(ix1-1:ix2+1,jy1-2:jy2,kz1) )/2./dy)/ &
        !     sqrt( ((s(ix1:ix2+2,jy1-1:jy2+1,kz1) - &
        !             s(ix1-2:ix2,jy1-1:jy2+1,kz1))/2./dx)**2 &
        !         + ((s(ix1-1:ix2+1,jy1:jy2+2,kz1) - &
        !             s(ix1-1:ix2+1,jy1-2:jy2,kz1))/2./dy)**2 )


        pf(ix1-1:ix2+1,jy1-1:jy2+1,k)   = 0.0
        pf(ix1-1:ix2+1,jy1-1:jy2+1,k)   = (sign(1.0,s(ix1-1:ix2+1,jy1-1:jy2+1,k))+1.0)/2.0


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

        subroutine mph_KPDcurvature2DC(s,crv,rho1x,rho2x,rho1y,rho2y,pf,w,sigx,sigy,dx,dy, &
           rho1,rho2,xit,crmx,crmn,ix1,ix2,jy1,jy2,thco1,thco2,cp1,cp2,mdot,tmic,lambda,blockID)   

   
        use Multiphase_data, ONLY : mph_meshMe
        use Heat_AD_data, only : ht_Tsat,ht_fmic

        use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_getBlkCenterCoords, Grid_getDeltas

        implicit none

#include "Flash.h"
#include "constants.h"

        !implicit real*8(a-h,o-z)
        !common/param/re,xit,rho1,rho2,g,sp,ubc,uout,nint

        integer, intent(in) :: ix1,ix2,jy1,jy2
        real, intent(in) :: dx, dy, rho1, rho2, xit ,thco1,thco2,cp1,cp2
        real, intent(out) :: crmx, crmn

        real, dimension(:,:,:), intent(inout):: s,crv,rho1x,rho2x,rho1y, &
                                               rho2y,pf,w,sigx,sigy,tmic

        real, dimension(:,:,:), intent(in) :: mdot,lambda

        integer, intent(in) :: blockID

        integer :: icrv(NXB+2*NGUARD,NYB+2*NGUARD,1)

        real :: ycell
        !- kpd - 
        real :: th,aa,xijl,xijr, &
                a1,a2,cri,xid,xij,xidl,xidr,yid,yidr,yidl,yij,yijl,yijr,bb,mfl,mfr,mT
        integer :: i,j,k
        real, parameter :: eps = 1E-13

        integer :: iSmear

        real del(MDIM),bsize(MDIM),coord(MDIM)
        real, dimension(2,MDIM) :: boundBox

        real :: pfl(NXB+2*NGUARD,NYB+2*NGUARD,1)
        real :: rho3, rhof

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
        w = 0.
        icrv = 0
        bb   = (rho2/rho1) - 1.0

        tmic = 0.

        !- kpd - Need to loop through one guard cell on each side to set jumps 
        !           when they cross block boundaries
                    !do j = jy1,jy2
                    !   do i = ix1,ix2
        k=1

        pfl(ix1-1:ix2+1,jy1-1:jy2+1,k)  = 0.0
        pfl(ix1-1:ix2+1,jy1-1:jy2+1,k)  = (sign(1.0,lambda(ix1-1:ix2+1,jy1-1:jy2+1,k))+1.0)/2.0

        rho3 = (rho1 + rho2)/2.0

        call Grid_getDeltas(blockID,del)
        call Grid_getBlkCenterCoords(blockId,coord)
        call Grid_getBlkBoundBox(blockId,boundBox)

        bsize(:) = boundBox(2,:) - boundBox(1,:)

        do j = jy1-1,jy2
           do i = ix1-1,ix2

                if(pfl(i,j,k) .eq. 0.0 .and. pfl(i+1,j,k) .eq. 1.0) then
                    th = abs(lambda(i+1,j,k))/(abs(lambda(i+1,j,k))+abs(lambda(i,j,k)))

                    if(pf(i,j,k) .eq. 0.0) then
                        rhof = rho2
                    else
                        rhof = rho1
                    end if

                    aa = th*(rho3/rho2) + (1-th)*(rhof/rho2)

                    rho2x(i+1,j,k) = rho2x(i+1,j,k)*(rho3/rho2)/aa
                end if
                !
                !
                if(pfl(i,j,k) .eq. 1.0 .and. pfl(i+1,j,k) .eq. 0.0) then
                    th = abs(lambda(i,j,k))/(abs(lambda(i+1,j,k))+abs(lambda(i,j,k)))

                    if(pf(i+1,j,k) .eq. 0.0) then
                        rhof = rho2
                    else
                        rhof = rho1
                    end if

                    aa = th*(rho3/rho2) + (1-th)*(rhof/rho2)

                    rho2x(i+1,j,k) = rho2x(i+1,j,k)*(rho3/rho2)/aa
                end if
                !
                !
                if(pfl(i,j,k) .eq. 0.0 .and. pfl(i,j+1,k) .eq. 1.0) then
                    th = abs(lambda(i,j+1,k))/(abs(lambda(i,j+1,k))+abs(lambda(i,j,k)))

                    if(pf(i,j,k) .eq. 0.0) then
                        rhof = rho2
                    else
                        rhof = rho1
                    end if

                    aa = th*(rho3/rho2) + (1-th)*(rhof/rho2)

                    rho2y(i,j+1,k) = rho2y(i,j+1,k)*(rho3/rho2)/aa
                end if
                !
                !
                if(pfl(i,j,k) .eq. 1.0 .and. pfl(i,j+1,k) .eq. 0.0) then
                    th = abs(lambda(i,j,k))/(abs(lambda(i,j+1,k))+abs(lambda(i,j,k)))

                    if(pf(i,j+1,k) .eq. 0.0) then
                        rhof = rho2
                    else
                        rhof = rho1
                    end if

                    aa = th*(rho3/rho2) + (1-th)*(rhof/rho2)

                    rho2y(i,j+1,k) = rho2y(i,j+1,k)*(rho3/rho2)/aa
                end if 

           end do
        end do

        do j = jy1-1,jy2
           do i = ix1-1,ix2


              ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                      real(j - NGUARD - 1)*del(JAXIS)  +  &
                      0.5*del(JAXIS)

              !--------------------------------------------------------------
              !- kpd - pf=0 (water) in current cell and pf=1 (air) in cell to right
              !--------------------------------------------------------------

              if(lambda(i,j,k) .lt. 0.0 .and. lambda(i+1,j,k) .lt. 0.0) then
              if(pf(i,j,k).eq.0..and.pf(i+1,j,k).eq.1. ) then

                 !          = (+)            = (+)           = (-)
                 th = abs(s(i+1,j,k))/(abs(s(i+1,j,k))+abs(s(i,j,k)))
                 
                 cri = crv(i+1,j,k)*(1.-th) + crv(i,j,k)*th
    
                 xijl = xit*crv(i,j,k)                 !- kpd - sigma*K. Used for jump in pressure
                 xijr = xit*crv(i+1,j,k)               !- kpd - sigma*K. Used for jump in pressure
                 xidl = 0.                             !- kpd - Used for jump in gradient
                 xidr = 0.                             !- kpd - Used for jump in gradient
                 mfl  = bb*mdot(i,j,k)*mdot(i,j,k)
                 mfr  = bb*mdot(i+1,j,k)*mdot(i+1,j,k)

                 xij = xijl*th + xijr*(1.-th)          !- kpd - Jump in value
                 xid = xidl*th + xidr*(1.-th)          !- kpd - Jump in gradient. Equal to 0 here.
                 mT  = mfl*th  + mfr*(1.-th)
                 
                 if (iSmear  .eq. 1) then

                 !- kpd - All Densities are relative to rho2...

                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)           !- kpd - Mixture density (not inverse)
                 rho1x(i+1,j,k) = rho1x(i+1,j,k)*(rho1/rho2)/aa
                 rho2x(i+1,j,k) = rho2x(i+1,j,k)*(rho2/rho2)/aa

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k)   - xij/aa/dx**2 - xid*th*(rho1/rho2)/aa/dx      + mT/aa/dx**2
                 w(i+1,j,k) = w(i+1,j,k) + xij/aa/dx**2 - xid*(1.-th)*(rho2/rho2)/aa/dx - mT/aa/dx**2
                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigx(i+1,j,k) = - xij/aa/dx + mT/aa/dx           !- kpd - sigma*K/rho/dx 
         
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
              if(pf(i,j,k).eq.1..and.pf(i+1,j,k).eq.0. ) then

                 th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k)))

                 cri = crv(i,j,k)*(1.-th) + crv(i+1,j,k)*th

                 xijl = xit*crv(i,j,k)
                 xijr = xit*crv(i+1,j,k)
                 xidl = 0. 
                 xidr = 0.
                 mfl  = bb*mdot(i,j,k)*mdot(i,j,k)
                 mfr  = bb*mdot(i+1,j,k)*mdot(i+1,j,k)

                 xij = xijl*(1.-th) + xijr*th
                 xid = xidl*(1.-th) + xidr*th
                 mT  = mfl*(1.-th)  + mfr*th

                 if (iSmear  .eq. 1) then

                 !- kpd - All Densities are relative to rho2...

                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)
                 rho1x(i+1,j,k) = rho1x(i+1,j,k)*(rho1/rho2)/aa
                 rho2x(i+1,j,k) = rho2x(i+1,j,k)*(rho2/rho2)/aa   

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k)   + xij/aa/dx**2 + xid*(1.-th)*(rho2/rho2)/aa/dx - mT/aa/dx**2
                 w(i+1,j,k) = w(i+1,j,k) - xij/aa/dx**2 + xid*th*(rho1/rho2)/aa/dx      + mT/aa/dx**2
                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigx(i+1,j,k) = xij/aa/dx - mT/aa/dx

                 else

                 end if

                 crmx = max(abs(cri),crmx)
                 crmn = min(abs(cri),crmn)

                 icrv(i,j,k) = 1
                 icrv(i+1,j,k) = 1

              end if

              end if

              if(lambda(i,j,k) .lt. 0.0 .and. lambda(i,j+1,k) .lt. 0.0) then

              !--------------------------------------------------------------
              !- kpd - pf=0 in current cell and pf=1 in cell above
              !--------------------------------------------------------------
              if(pf(i,j,k).eq.0..and.pf(i,j+1,k).eq.1. ) then

                 th = abs(s(i,j+1,k))/(abs(s(i,j+1,k))+abs(s(i,j,k)))

                 cri = crv(i,j+1,k)*(1.-th) + crv(i,j,k)*th

                 yijl = xit*crv(i,j,k)
                 yijr = xit*crv(i,j+1,k)
                 yidl = 0. 
                 yidr = 0.
                 mfl  = bb*mdot(i,j,k)*mdot(i,j,k)
                 mfr  = bb*mdot(i,j+1,k)*mdot(i,j+1,k)

                 yij = yijl*th + yijr*(1.-th)
                 yid = yidl*th + yidr*(1.-th)
                 mT  = mfl*th  + mfr*(1.-th)

                 if (iSmear  .eq. 1) then

                 !- kpd - All Densities are relative to rho2...

                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)
                 rho1y(i,j+1,k) = rho1y(i,j+1,k)*(rho1/rho2)/aa
                 rho2y(i,j+1,k) = rho2y(i,j+1,k)*(rho2/rho2)/aa 

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k) - yij/aa/dy**2 - yid*th*(rho1/rho2)/aa/dy         + mT/aa/dy**2
                 w(i,j+1,k) = w(i,j+1,k)   + yij/aa/dy**2 - yid*(1.-th)*(rho2/rho2)/aa/dy- mT/aa/dy**2
                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigy(i,j+1,k) = - yij/aa/dy + mT/aa/dy

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
              if(pf(i,j,k).eq.1..and.pf(i,j+1,k).eq.0. ) then

                 th = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k))) 
 
                 cri = crv(i,j,k)*(1.-th) + crv(i,j+1,k)*th

                 yijl = xit*crv(i,j,k)
                 yijr = xit*crv(i,j+1,k)
                 yidl = 0. 
                 yidr = 0.
                 mfl  = bb*mdot(i,j,k)*mdot(i,j,k)
                 mfr  = bb*mdot(i,j+1,k)*mdot(i,j+1,k)

                 yij = yijl*(1.-th) + yijr*th
                 yid = yidl*(1.-th) + yidr*th
                 mT  = mfl*(1.-th)  + mfr*th

                 if (iSmear  .eq. 1) then

                 !- kpd - All Densities are relative to rho2...

                 aa = th*(rho1/rho2) + (1.-th)*(rho2/rho2)
                 rho1y(i,j+1,k) = rho1y(i,j+1,k)*(rho1/rho2)/aa
                 rho2y(i,j+1,k) = rho2y(i,j+1,k)*(rho2/rho2)/aa  

                 !- kpd - "w" is the Source term for Pressure Eqn
                 w(i,j,k)   = w(i,j,k)   + yij/aa/dy**2 + yid*(1.-th)*(rho2/rho2)/aa/dy - mT/aa/dy**2
                 w(i,j+1,k) = w(i,j+1,k) - yij/aa/dy**2 + yid*th*(rho1/rho2)/aa/dy      + mT/aa/dy**2
                 !- kpd - "sig" is the source term in Momentum Equations. Only uses 
                 !           the jump in value, not the jump in derivative.
                 sigy(i,j+1,k) = yij/aa/dy - mT/aa/dy

                 else

                 end if

                 crmx = max(abs(cri),crmx)
                 crmn = min(abs(cri),crmn)

                 icrv(i,j,k) = 1
                 icrv(i,j+1,k) = 1

              end if
              end if

           end do
        end do

        !- kpd - This is done for post-processing to visualize where jumps were applied...
        crv = crv*icrv

      end subroutine mph_KPDcurvature2DC
