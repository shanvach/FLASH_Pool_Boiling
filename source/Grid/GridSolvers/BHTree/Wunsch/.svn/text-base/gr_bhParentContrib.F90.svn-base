!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhParentContrib
!!
!! NAME
!!
!!  gr_bhParentContrib
!!
!!
!! SYNOPSIS
!!
!!   gr_bhParentContrib(
!!          integer,intent(in) :: block,
!!          integer,intent(in) :: tr,
!!          integer,intent(in) :: cpu,
!!          real,intent(INOUT),dimension(1:gr_bhTreeBS, 1:gr_bhTreeBS, 1:gr_bhTreeBS) :: Phi
!!          )
!!
!! DESCRIPTION
!!
!!   Computes the contribution of a parent block to the gravitational potential at a specific point.
!!
!! ARGUMENTS
!!
!!  block   - ID of a block where the potential is calculated
!!  tr      - ID of a block/tree which contributes to the potential
!!  cpu     - cpu of a block/tree which contributes to the potential
!!  Phi     - 3D array containing potention, contribution added to it
!!
!!***

subroutine gr_bhParentContrib(block, tr, cpu, Phi)

  use gr_bhInterface, ONLY : gr_bhEwald, gr_bhGetTreeSize
  use gr_bhData, ONLY : GR_TREE_BND_PERIODIC, gr_bhBndType, &
    gr_bhLx, gr_bhLy, gr_bhLz, gr_bhTreeBS, gr_bhLocCoords, &
    gr_bhTreeParentTree

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer,intent(in) :: block, tr, cpu
  real, intent(INOUT), dimension(1:gr_bhTreeBS, 1:gr_bhTreeBS, 1:gr_bhTreeBS) :: Phi
  integer :: ii, jj, kk
  real  ::  disti, dx, dy, dz

  ! block is far enough to be added to the potential
  do kk = 1,gr_bhTreeBS ! this block
    do jj = 1,gr_bhTreeBS
      do ii = 1,gr_bhTreeBS
 
        ! determine the distance between this cell and the parent block (tr, cpu)
        dx = (gr_bhLocCoords(ii, 1, block) - gr_bhTreeParentTree(2, tr, cpu))
        dy = (gr_bhLocCoords(jj, 2, block) - gr_bhTreeParentTree(3, tr, cpu))
        dz = (gr_bhLocCoords(kk, 3, block) - gr_bhTreeParentTree(4, tr, cpu))

        ! find the minimum distance with periodic boundaries
        if (gr_bhBndType(1) .EQ. GR_TREE_BND_PERIODIC) then
          dx = min(abs(dx), abs(dx+gr_bhLx), abs(dx-gr_bhLx))
        endif
        if (gr_bhBndType(3) .EQ. GR_TREE_BND_PERIODIC) then
          dy = min(abs(dy), abs(dy+gr_bhLy), abs(dy-gr_bhLy))
        endif
        if (gr_bhBndType(5) .EQ. GR_TREE_BND_PERIODIC) then
          dz = min(abs(dz), abs(dz+gr_bhLz), abs(dz-gr_bhLz))
        endif
          
        disti = 1.0/sqrt(dx*dx + dy*dy + dz*dz + 1.d-99)
          
        if ( (gr_bhBndType(1) .EQ. GR_TREE_BND_PERIODIC) &
        .or. (gr_bhBndType(3) .EQ. GR_TREE_BND_PERIODIC) &
        .or. (gr_bhBndType(5) .EQ. GR_TREE_BND_PERIODIC)) then
          Phi(ii,jj,kk) = Phi(ii,jj,kk) &
          & + gr_bhTreeParentTree(1, tr, cpu) * gr_bhEwald(abs(dx), abs(dy), abs(dz))
        else
          Phi(ii,jj,kk) = Phi(ii,jj,kk) &
          & + gr_bhTreeParentTree(1, tr, cpu) * disti
        endif
      enddo
    enddo
  enddo

end subroutine gr_bhParentContrib
