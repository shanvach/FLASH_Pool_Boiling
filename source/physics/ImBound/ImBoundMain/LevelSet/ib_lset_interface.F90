Module ib_lset_interface

  implicit none

#include "constants.h"
#include "Flash.h"
#include "IncompNS.h"

interface
      subroutine ib_lset(blockCount,blockList,dt)

      integer, INTENT(INOUT) :: blockCount
      integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList
      real,    INTENT(IN) :: dt

      end subroutine
end interface

interface
     subroutine ib_advectWENO3(s,u,v,dt,dx,dy,ix1,ix2,jy1,jy2)
        real, dimension(:,:,:), intent(inout):: s
        real, dimension(:,:,:), intent(in) :: u,v
        real, intent(in) :: dt,dx,dy
        integer, intent(in) :: ix1,ix2,jy1,jy2
     end subroutine
end interface

interface
        subroutine ib_advectWENO3_3D(s,u,v,w,dt,dx,dy,dz,&
                                        ix1,ix2,jy1,jy2,kz1,kz2)
        real, dimension(:,:,:), intent(inout):: s
        real, dimension(:,:,:), intent(in)   :: u,v,w
        real, intent(in) :: dt,dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        end subroutine
end interface

interface
     subroutine ib_advect(blockCount,blockList,dt)
      integer, INTENT(INOUT) :: blockCount
      integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList
      real,    INTENT(IN) :: dt
     end subroutine
end interface

interface
      subroutine ib_lset_3D(blockCount,blockList,dt)
      integer, INTENT(INOUT) :: blockCount
      integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList
      real,    INTENT(IN) :: dt
      end subroutine
end interface

interface
      subroutine ib_lsRedistance(s,u,v,dx,dy,ix1,ix2,jy1,jy2,soo,lsDT, blockID,minCellDiag)
        integer, intent(in) :: ix1,ix2,jy1,jy2,blockID
        real, dimension(:,:,:), intent(inout):: s
        real, intent(in) :: dx,dy,lsDT, minCellDiag
        real, dimension(:,:,:), intent(in):: u,v,soo
      end subroutine
end interface

interface
      subroutine ib_lsRedistance_3D(s,u,v,w,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,soo,lsDT,minCellDiag)
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
        real, dimension(:,:,:), intent(inout):: s
        real, intent(in) :: dx,dy,dz,lsDT, minCellDiag
        real, dimension(:,:,:), intent(in):: u,v,w,soo
      end subroutine
end interface

End module ib_lset_interface
