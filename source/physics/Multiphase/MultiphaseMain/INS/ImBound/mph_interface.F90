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
        subroutine mph_KPDcurvature2DAB(s,lambda,crv,rho1x,rho2x,rho1y,rho2y,pf,w,sigx,sigy,dx,dy, &
           rho1,rho2,xit,crmx,crmn,ix1,ix2,jy1,jy2,visc,vis1,vis2,blockID)
        implicit none
        integer, intent(in) :: ix1,ix2,jy1,jy2,blockID
        real, intent(in) :: dx, dy, rho1, rho2, xit, vis1, vis2
        real, intent(out) :: crmx, crmn

        real, dimension(:,:,:), intent(inout):: s,crv,rho1x,rho2x,rho1y, &
                                                rho2y,pf,w,sigx,sigy,visc,&
                                                lambda
        end subroutine mph_KPDcurvature2DAB
end interface

interface
        subroutine mph_KPDcurvature2DC(s,lambda,crv,rho1x,rho2x,rho1y,rho2y,pf,w,sigx,sigy,dx,dy, &
           rho1,rho2,xit,crmx,crmn,ix1,ix2,jy1,jy2,blockID)
        implicit none
        integer, intent(in) :: ix1,ix2,jy1,jy2
        real, intent(in) :: dx, dy, rho1, rho2, xit
        real, intent(out) :: crmx, crmn

        real, dimension(:,:,:), intent(inout):: s,crv,rho1x,rho2x,rho1y, &
                                                rho2y,pf,w,sigx,sigy,lambda
        integer, intent(in) :: blockID
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
        subroutine mph_KPDcurvature3DAB(s,lambda,crv,dx,dy,dz, &
           ix1,ix2,jy1,jy2,kz1,kz2, &
           rho1x,rho2x,rho1y,rho2y,rho1z,rho2z,pf,rho1,rho2,visc,vis1,vis2,blockID)

        implicit none

        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2,blockID
        real, intent(in) :: dx, dy, dz, rho1, rho2, vis1, vis2
        real, dimension(:,:,:), intent(inout):: s,crv, &
                                                rho1x,rho2x,rho1y, &
                                                rho2y,pf, &
                                                rho1z,rho2z, visc,lambda

        end subroutine mph_KPDcurvature3DAB
end interface

interface
        subroutine mph_KPDcurvature3DC(s,lambda,crv,rho1x,rho2x,rho1y,rho2y, &
                                       pf,w,sigx,sigy,dx,dy,          &
                                       rho1,rho2,xit,ix1,ix2, &
                                       jy1,jy2,dz,kz1,kz2,rho1z, &
                                       rho2z,sigz,blockID)
        implicit none
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        real, intent(in) :: dx, dy, dz, rho1, rho2, xit

        real, dimension(:,:,:), intent(inout):: s,crv,rho1x,rho2x,rho1y, &
                                                rho2y,pf,w,sigx,sigy, &
                                                rho1z,rho2z,sigz,lambda

        integer, intent(in) :: blockID
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
      integer, INTENT(IN) :: mph_flag 
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
      subroutine mph_imbound(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

      integer, INTENT(IN) :: sweepOrder
      integer, INTENT(INOUT) :: blockCount
      integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList !blockCount
      real,    INTENT(IN) :: timeEndAdv,dt,dtOld

      end subroutine
end interface

interface
      subroutine mph_imboundExtrap(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

      integer, INTENT(IN) :: sweepOrder
      integer, INTENT(INOUT) :: blockCount
      integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList !blockCount
      real,    INTENT(IN) :: timeEndAdv,dt,dtOld

      end subroutine
end interface
! End

End module mph_interface
