!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhILContrib
!!
!! NAME
!!
!!  gr_bhILContrib
!!
!!
!! SYNOPSIS
!!
!!   gr_bhILContrib(
!!          integer,intent(in) :: block, tr, cpu, temp
!!          real,intent(INOUT),dimension(1:gr_bhTreeBS, 1:gr_bhTreeBS, 1:gr_bhTreeBS) :: Phi
!!          real,intent(out)   :: cellcnt, nodecnt
!!          )
!!
!! DESCRIPTION
!!
!!   Computes the contribution of a leaf block to the gravitational potential at a specific point.
!!   Uses interaction lists.
!!
!! ARGUMENTS
!!  
!!  block   - ID of a block where the potential is calculated
!!  tr      - ID of a block/tree which contributes to the potential
!!  cpu     - cpu of a block/tree which contributes to the potential
!!  temp    - template number (2-27; one of 26 3D-neighbours)
!!  Phi     - 3D array containing potention, contribution added to it
!!  cellcnt - number of cell-cell interactions
!!  nodecnt - number of cell-node interactions
!!
!!***

subroutine gr_bhILContrib(block, tr, cpu, temp, Phi, cellcnt, nodecnt)

  use gr_bhInterface, ONLY : gr_bhEwald, gr_bhGetTreeSize
  use gr_bhData, ONLY : GR_TREE_BND_PERIODIC, gr_bhBndType, &
    gr_bhLx, gr_bhLy, gr_bhLz, gr_bhTreeArray, gr_bhTreeBS, &
    gr_bhLocCoords, gr_bhTreeMyPE, gr_bhTreeLrefine, gr_bhTreeCellSize, &
    gr_bhTreeCc_count, gr_bhTreeCc_disti, gr_bhTreeCn_ind, gr_bhTreeCc_ind

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer,intent(in) :: block, tr, cpu, temp
  real,intent(INOUT),dimension(1:gr_bhTreeBS, 1:gr_bhTreeBS, 1:gr_bhTreeBS) :: Phi
  real,intent(out)   :: cellcnt, nodecnt
  integer :: ii, jj, kk
  integer :: i, j, k, cind
  real  :: mass, disti, dx, dy, dz, dx_celli
  real  :: x, y, z, xx, yy, zz

  cellcnt = 0.0 ! cell counter
  nodecnt = 0.0 ! node counter


  ! get size of a cell in x-direction
  dx_celli = 1./gr_bhTreeCellSize(gr_bhTreeLrefine(block, gr_bhTreeMyPE), IAXIS)

  ! go through all cells in the block
  do ii = 1,gr_bhTreeBS
    do jj = 1,gr_bhTreeBS
      do kk = 1,gr_bhTreeBS

        xx = gr_bhLocCoords(ii, 1, block)
        yy = gr_bhLocCoords(jj, 2, block)
        zz = gr_bhLocCoords(kk, 3, block)
        cind = (ii-1) + (jj-1)*gr_bhTreeBS + (kk-1)*gr_bhTreeBS**2 + 1
        ! add cells
        i = 1
        do
          if (gr_bhTreeCc_count(cind,i,temp) .eq. -1) exit ! check for the end of the list

          ! sum the mass of cells at the same distance
          mass = 0
          do j = 1, gr_bhTreeCc_count(cind,i,temp)
            mass = mass + gr_bhTreeArray(cpu, tr)%p(gr_bhTreeCc_ind(cind,i,j,temp))
          enddo

          ! read the distance and add to the potential
          disti = gr_bhTreeCc_disti(cind,i,temp)*dx_celli
          Phi(ii,jj,kk) = Phi(ii,jj,kk) + mass * disti

          i = i + 1
          cellcnt = cellcnt + 1
        enddo

        ! add nodes
        i = 1
        do
          if (gr_bhTreeCn_ind(cind,i,temp) .eq. -1) exit

          ! read the mass centre coords from the gr_bhTreeArray
          x = gr_bhTreeArray(cpu, tr)%p(gr_bhTreeCn_ind(cind,i,temp)+1)
          y = gr_bhTreeArray(cpu, tr)%p(gr_bhTreeCn_ind(cind,i,temp)+2)
          z = gr_bhTreeArray(cpu, tr)%p(gr_bhTreeCn_ind(cind,i,temp)+3)

          ! compute square of distance
          dx = x - xx
          dy = y - yy
          dz = z - zz
          if (gr_bhBndType(1) .EQ. GR_TREE_BND_PERIODIC) then
            dx = min(abs(dx), abs(dx+gr_bhLx), abs(dx-gr_bhLx))
          endif
          if (gr_bhBndType(3) .EQ. GR_TREE_BND_PERIODIC) then
            dy = min(abs(dy), abs(dy+gr_bhLy), abs(dy-gr_bhLy))
          endif
          if (gr_bhBndType(5) .EQ. GR_TREE_BND_PERIODIC) then
            dz = min(abs(dz), abs(dz+gr_bhLz), abs(dz-gr_bhLz))
          endif
          disti = sqrt(1.0 / (dx*dx + dy*dy + dz*dz + 1.d-99))

          Phi(ii,jj,kk) = Phi(ii,jj,kk) &
          & + gr_bhTreeArray(cpu, tr)%p(gr_bhTreeCn_ind(cind,i,temp)) * disti

          i = i + 1
          nodecnt = nodecnt + 1
        enddo

      enddo
    enddo
  enddo

end subroutine gr_bhILContrib

