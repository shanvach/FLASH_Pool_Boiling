!!****if* source/physics/Multiphase/MultiphaseMain/mph_interface
!!
!! NAME
!!  
!!
!! SYNOPSIS
!!  mph_interface()
!!
!! DESCRIPTION
!!  This is an interface specific for the Multiphase Incompressible Navier Stokes
!!  module that defines its public interfaces.
!!
!!***
Module mph_interface

  implicit none

#include "constants.h"
#include "Flash.h"
#include "IncompNS.h"


!interface
!        subroutine mph_KPDcurvature2DA(s,crv,rho1x,rho2x,rho1y,rho2y,pf,w,sigx,sigy,dx,dy, &
!           rho1,rho2,xit,crmx,crmn,ix1,ix2,jy1,jy2)   
!        implicit none
!        integer, intent(in) :: ix1,ix2,jy1,jy2
!        real, intent(in) :: dx, dy, rho1, rho2, xit 
!        real, intent(out) :: crmx, crmn
!
!        real, dimension(:,:,:), intent(inout):: s,crv,rho1x,rho2x,rho1y, &
!                                                rho2y,pf,w,sigx,sigy
!        end subroutine mph_KPDcurvature2DA
!end interface
!
!interface
!        subroutine mph_KPDcurvature2DB(s,crv,rho1x,rho2x,rho1y,rho2y,pf,w,sigx,sigy,dx,dy, &
!           rho1,rho2,xit,crmx,crmn,ix1,ix2,jy1,jy2,visc)
!        implicit none
!        integer, intent(in) :: ix1,ix2,jy1,jy2
!        real, intent(in) :: dx, dy, rho1, rho2, xit
!        real, intent(out) :: crmx, crmn
!
!        real, dimension(:,:,:), intent(inout):: s,crv,rho1x,rho2x,rho1y, &
!                                                rho2y,pf,w,sigx,sigy,visc
!        end subroutine mph_KPDcurvature2DB
!end interface

interface
        subroutine mph_KPDcurvature2DAB(s,crv,rho1x,rho2x,rho1y,rho2y,pf,w,sigx,sigy,dx,dy, &
           rho1,rho2,xit,crmx,crmn,ix1,ix2,jy1,jy2,visc,vis1,vis2,alph,thco1,thco2,cp1,cp2,nrmx,nrmy,mflg,smhv,smrh,&
           lambda)
        implicit none
        integer, intent(in) :: ix1,ix2,jy1,jy2
        real, intent(in) :: dx,dy,rho1,rho2,xit,vis1,vis2,thco1,thco2,cp1,cp2
        real, intent(out) :: crmx, crmn

        real, dimension(:,:,:), intent(inout):: lambda,s,crv,rho1x,rho2x,rho1y, &
                                                rho2y,pf,w,sigx,sigy,visc,alph,nrmx,nrmy,mflg,smhv,smrh
        end subroutine mph_KPDcurvature2DAB
end interface

interface
        subroutine mph_KPDcurvature2DC(s,crv,rho1x,rho2x,rho1y,rho2y,pf,w,sigx,sigy,dx,dy, &
           rho1,rho2,xit,crmx,crmn,ix1,ix2,jy1,jy2,thco1,thco2,cp1,cp2,mdot,tmic,lambda,blockID)
        implicit none
        integer, intent(in) :: ix1,ix2,jy1,jy2,blockID
        real, intent(in) :: dx, dy, rho1, rho2, xit,thco1,thco2,cp1,cp2
        real, intent(out) :: crmx, crmn

        real, dimension(:,:,:), intent(inout):: s,crv,rho1x,rho2x,rho1y, &
                                                rho2y,pf,w,sigx,sigy,tmic

        real, dimension(:,:,:), intent(in) :: mdot,lambda
        end subroutine mph_KPDcurvature2DC
end interface

!interface
!        subroutine mph_KPDcurvature3DA(s,crv,dx,dy, &
!           ix1,ix2,jy1,jy2,dz,kz1,kz2)
!        implicit none
!        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
!        real, intent(in) :: dx, dy, dz
!
!        real, dimension(:,:,:), intent(inout):: s,crv
!
!        end subroutine mph_KPDcurvature3DA
!end interface
!
!interface
!        subroutine mph_KPDcurvature3DB(s,rho1x,rho2x,rho1y,rho2y,pf, &
!           rho1,rho2,ix1,ix2,jy1,jy2,kz1,kz2,rho1z,rho2z)
!
!        implicit none
!        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
!        real, intent(in) :: rho1, rho2                 
!
!        real, dimension(:,:,:), intent(inout):: s,rho1x,rho2x,rho1y, &
!                                                rho2y,pf, &
!                                                rho1z,rho2z
!        end subroutine mph_KPDcurvature3DB
!end interface

interface
        subroutine mph_KPDcurvature3DAB(s,crv,dx,dy,dz, &
           ix1,ix2,jy1,jy2,kz1,kz2, &
           rho1x,rho2x,rho1y,rho2y,rho1z,rho2z,pf,rho1,rho2,visc,vis1,vis2,alph,&
           thco1,thco2,cp1,cp2,nrmx,nrmy,nrmz,mflg,smhv,lambda)

        implicit none

        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        real, intent(in) :: dx, dy, dz, rho1, rho2, vis1, vis2,thco1,thco2,cp1,cp2
        real, dimension(:,:,:), intent(inout):: s,crv, &
                                                rho1x,rho2x,rho1y, &
                                                rho2y,pf, &
                                                rho1z,rho2z, visc,alph,&
                                                nrmx,nrmy,nrmz,mflg,smhv
        real, dimension(:,:,:), intent(in) :: lambda
        end subroutine mph_KPDcurvature3DAB
end interface

interface
        subroutine mph_KPDcurvature3DC(s,crv,rho1x,rho2x,rho1y,rho2y, &
                                       pf,w,sigx,sigy,dx,dy,          &
                                       rho1,rho2,xit,ix1,ix2, &
                                       jy1,jy2,dz,kz1,kz2,rho1z, &
                                       rho2z,sigz,mdot,tmic,temp,lambda,blockID)
        implicit none
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2,blockID
        real, intent(in) :: dx, dy, dz, rho1, rho2, xit

        real, dimension(:,:,:), intent(inout):: s,crv,rho1x,rho2x,rho1y, &
                                                rho2y,pf,w,sigx,sigy, &
                                                rho1z,rho2z,sigz,tmic

        real, dimension(:,:,:), intent(in) :: mdot,temp,lambda

        end subroutine mph_KPDcurvature3DC
end interface

interface
        subroutine mph_KPDadvectWENO3(s,u,v,dt,dx,dy,ix1,jy1,ix2,jy2,lambda,blockID) 

        real, dimension(:,:,:), intent(inout):: s,u,v,lambda
        real, intent(in) :: dt,dx,dy
        integer, intent(in) :: ix1,ix2,jy1,jy2,blockID

        end subroutine mph_KPDadvectWENO3
end interface

interface
        subroutine mph_KPDadvectWENO3_3D(s,u,v,w,dt,dx,dy,dz,ix1,jy1,ix2,jy2,kz1,kz2,lambda,blockID)

        real, dimension(:,:,:), intent(inout):: s,u,v,w,lambda
        real, intent(in) :: dt,dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2,blockID

        end subroutine mph_KPDadvectWENO3_3D
end interface


interface
        subroutine mph_KPDlsRedistance(s,u,v,dx,dy,ix1,ix2,jy1,jy2,soo,lsDT,blockID,minCellDiag)

        integer, intent(in) :: ix1,ix2,jy1,jy2,blockID
        real, dimension(:,:,:), intent(in):: u,v,soo
        real, dimension(:,:,:), intent(inout):: s
        real, intent(in) :: dx,dy,lsDT,minCellDiag

        end subroutine mph_KPDlsRedistance
end interface

interface
        subroutine mph_KPDlsRedistance_3D(s,u,v,w,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,soo,lsDT,minCellDiag)

        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        real, dimension(:,:,:), intent(in):: u,v,w,soo
        real, dimension(:,:,:), intent(inout):: s
        real, intent(in) :: dx,dy,dz,lsDT, minCellDiag

        end subroutine mph_KPDlsRedistance_3D
end interface

! Start
! This subroutine is added by Akash
interface
      subroutine mph_evolve(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder,mph_flag)

      integer, INTENT(IN) :: sweepOrder
      integer, INTENT(INOUT) :: blockCount
      integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList !blockCount
      real,    INTENT(IN) :: timeEndAdv,dt,dtOld
      integer, intent(in) :: mph_flag
       
      end subroutine
end interface

interface
      subroutine mph_advect(blockCount, blockList,timeEndAdv,dt,dtOld,sweepOrder)

      integer, INTENT(IN) :: sweepOrder
      integer, INTENT(INOUT) :: blockCount
      integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList !blockCount
      real,    INTENT(IN) :: timeEndAdv,dt,dtOld

      end subroutine
end interface

interface
        subroutine mph_getSmearedProperties2D(s,pf,dx,dy,rho1,rho2,ix1,ix2,jy1,jy2,xnorm,ynorm,smhv,smrh,lambda)

        integer, intent(in) :: ix1,ix2,jy1,jy2
        real, intent(in) :: dx,dy,rho1,rho2
        real, dimension(:,:,:), intent(in) :: xnorm,ynorm,s,pf,lambda
        real, dimension(:,:,:), intent(inout) :: smhv,smrh

        end subroutine
end interface

interface
       subroutine mph_getSmearedProperties3D(s,pf,dx,dy,dz,rho1,rho2,ix1,ix2,jy1,jy2,kz1,kz2,xnorm,ynorm,znorm,smhv,smrh,lambda) 

        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        real, intent(in) :: dx,dy,rho1,rho2,dz
        real, dimension(:,:,:), intent(in) :: xnorm,ynorm,s,znorm,pf,lambda
        real, dimension(:,:,:), intent(inout) :: smhv,smrh

      end subroutine
end interface


interface
        subroutine mph_getInterfaceVelocity(uni,vni,ru1,ix1,ix2,jy1,jy2,dx,dy, &
                           visc,rho1x,rho2x,rho1y,rho2y,gravX,gravY, &
                           mdot,smrh,xnorm,ynorm,uint,vint,curv)

        implicit none
        INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2
        REAL, INTENT(IN):: ru1, dx, dy, gravX,gravY
        REAL, DIMENSION(:,:,:), INTENT(IN):: uni, vni, visc, rho1x, rho2x, rho1y, rho2y
        REAL, DIMENSION(:,:,:), INTENT(IN) :: xnorm,ynorm,mdot,smrh,curv
        REAL, DIMENSION(:,:,:), INTENT(OUT):: uint, vint

        end subroutine
end interface

interface
        subroutine mph_getInterfaceVelocity_3D(uni,vni,wni,ru1,ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz,&
                           visc,rho1x,rho2x,rho1y,rho2y,rho1z,rho2z,gravX,gravY,gravZ,&
                           mdot,smrh,xnorm,ynorm,znorm,uint,vint,wint,curv)

        implicit none
        INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2,kz1,kz2
        REAL, INTENT(IN):: ru1, dx, dy, dz, gravX, gravY, gravZ
        REAL, DIMENSION(:,:,:), INTENT(IN):: uni, vni, wni, visc, rho1x, rho2x, rho1y, rho2y, rho1z, rho2z
        REAL, DIMENSION(:,:,:), INTENT(IN) :: xnorm,ynorm,znorm,mdot,smrh,curv
        REAL, DIMENSION(:,:,:), INTENT(OUT):: uint, vint, wint

        end subroutine
end interface

interface
       subroutine mph_extrapVars(u,v,T,p,s,pf,dx,dy,dz,nx,ny,ix1,ix2,jy1,jy2)
       implicit none
       real, dimension(:,:,:), intent(inout) :: u,v
       real, dimension(:,:,:), intent(in) :: T,p,s,pf,nx,ny
       real, intent(in) :: dx,dy,dz
       integer, intent(in) :: ix1,ix2,jy1,jy2
       end subroutine
end interface

interface
      subroutine mph_imbound(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

      integer, INTENT(IN) :: sweepOrder
      integer, INTENT(INOUT) :: blockCount
      integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList !blockCount
      real,    INTENT(IN) :: timeEndAdv,dt,dtOld

      end subroutine
end interface

interface
      subroutine mph_iblset(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)
      integer, INTENT(IN) :: sweepOrder
      integer, INTENT(INOUT) :: blockCount
      integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList
      real,    INTENT(IN) :: timeEndAdv,dt,dtOld
      end subroutine
end interface
! End

End module mph_interface
