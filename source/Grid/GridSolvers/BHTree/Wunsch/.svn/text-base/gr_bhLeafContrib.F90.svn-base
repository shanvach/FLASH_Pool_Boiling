!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhLeafContrib
!!
!! NAME
!!
!!  gr_bhLeafContrib
!!
!!
!! SYNOPSIS
!!
!!   gr_bhLeafContrib(
!!          integer,intent(in) :: block, tr, cpu, temp
!!          real,intent(INOUT),dimension(1:gr_bhTreeBS, 1:gr_bhTreeBS, 1:gr_bhTreeBS) :: Phi
!!          real,intent(out)   :: cellcnt, nodecnt
!!          )
!!
!! DESCRIPTION
!!
!!   Computes the contribution of a leaf block to the gravitational potential at a specific point.
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



subroutine gr_bhLeafContrib(block, tr, cpu, Phi, cellcnt, nodecnt)

  use gr_bhInterface, ONLY : gr_bhEwald, gr_bhGetTreeSize
  use gr_bhData, ONLY : GR_TREE_BND_PERIODIC, gr_bhBndType, &
    gr_bhLx, gr_bhLy, gr_bhLz, gr_bhTreeLevels, &
    gr_bhTreeArray, gr_bhTreeBS, gr_bhTreeLimAngle2, &
    gr_bhLocCoords, gr_bhTreeLoff, gr_bhTreeMyPE, &
    gr_bhTreeBCen, gr_bhTreeLrefine, gr_bhTreeCellSize, &
    gr_bhTreeDiag2

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"


  integer    :: ii, jj, kk
  integer,intent(in) :: block, tr, cpu
  real,intent(out)   :: cellcnt, nodecnt
  real, intent(INOUT), dimension(1:gr_bhTreeBS, 1:gr_bhTreeBS, 1:gr_bhTreeBS) :: Phi
  integer       :: i, j, k, l, level, sp
  integer       :: multi(1:gr_bhTreeLevels), mi(1:gr_bhTreeLevels)
  integer :: fs, fac, pos
  integer :: stack(1:gr_bhTreeLevels, 1:(8**(gr_bhTreeLevels+1)-1)/7)
  real  :: potential, mass, dist2i, dx, dy, dz
  real  :: x, y, z, xx, yy, zz


  cellcnt = 0.0 ! "distance to cell" counter
  nodecnt = 0.0 ! "distance to tree node" counter
  do kk = 1,gr_bhTreeBS ! this block
    do jj = 1,gr_bhTreeBS
      do ii = 1,gr_bhTreeBS

        potential = 0.0
   
        ! coords of the cell for which the potential is computed
        xx = gr_bhLocCoords(ii, 1, block)
        yy = gr_bhLocCoords(jj, 2, block)
        zz = gr_bhLocCoords(kk, 3, block)
   
        ! put root node on stack
        do l = 1, gr_bhTreeLevels
          multi(l) = 0
        enddo
        stack(:,1) = multi(:)
        sp = 1
   
        ! walk through the tree
        do
          ! take the node on the bottom of the stack
          multi(:) = stack(:,sp)
          sp = sp - 1
   
          ! determine level of the multi-index
          level = gr_bhTreeLevels
          do l = 1,gr_bhTreeLevels
            if (multi(l) == 0) then
              level = l - 1
              exit
            endif
          enddo
          
          ! get position of the node mass center
          if (level == gr_bhTreeLevels) then
            ! the lowest level, the node is just one cell
            ! at first convert multi-index to indeces in the block
            i = 1
            j = 1
            k = 1
            do l = 1,gr_bhTreeLevels
              fac = 2**(gr_bhTreeLevels-l)
              i = i + (mod(multi(l)-1,2))      * fac
              j = j + (mod((multi(l)-1)/2,2))  * fac
              k = k + (mod((multi(l)-1)/4,2))  * fac
            enddo
            
            if ((gr_bhTreeMyPE .eq. cpu) .and. (block .eq. tr)) then
              if ((i .eq. ii) .and. (j .eq. jj) .and. (k .eq. kk)) goto 10 ! omit the same cell
            endif
                
            ! and then calculate the coords from the global gr_bhTreeBCen array
            x = gr_bhTreeBCen(IAXIS, tr, cpu) - 0.5*gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), IAXIS)*(gr_bhTreeBS-2*i+1)
            y = gr_bhTreeBCen(JAXIS, tr, cpu) - 0.5*gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), JAXIS)*(gr_bhTreeBS-2*j+1)
            z = gr_bhTreeBCen(KAXIS, tr, cpu) - 0.5*gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), KAXIS)*(gr_bhTreeBS-2*k+1)
            !print *, 'xyzc = ', x, y, z
   
            ! compute position in the gr_bhTreeArray and get the mass
            fs = 1
            pos = gr_bhTreeLoff(level)
            do l = level,1,-1
              pos = pos + (multi(l)-1)*fs
              fs = fs * 8
            enddo
            mass = gr_bhTreeArray(cpu, tr)%p(pos)
   
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
            dist2i = 1.0 / (dx*dx + dy*dy + dz*dz + 1.d-99)
            cellcnt = cellcnt + 1 ! update the distance counter

            ! add the contribution to the potential
            if ( (gr_bhBndType(1) .EQ. GR_TREE_BND_PERIODIC) &
            .or. (gr_bhBndType(3) .EQ. GR_TREE_BND_PERIODIC) &
            .or. (gr_bhBndType(5) .EQ. GR_TREE_BND_PERIODIC)) then
              potential = potential + mass * gr_bhEwald(abs(dx), abs(dy), abs(dz))
            else
              potential = potential + mass*sqrt(dist2i)
            endif
            
          else
            ! find the position in the gr_bhTreeArray
            fs = 4
            pos = gr_bhTreeLoff(level)
            do l = level,1,-1
              pos = pos + (multi(l)-1)*fs
              fs = fs * 8
            enddo
            
            ! read the mass centre coords from the gr_bhTreeArray
            x = gr_bhTreeArray(cpu, tr)%p(pos+1)
            y = gr_bhTreeArray(cpu, tr)%p(pos+2)
            z = gr_bhTreeArray(cpu, tr)%p(pos+3)

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
            dist2i = 1.0 / (dx*dx + dy*dy + dz*dz + 1.d-99)
   
            ! test size/dist < angle
            if (gr_bhTreeDiag2(level+gr_bhTreeLrefine(tr,cpu))*dist2i .lt. gr_bhTreeLimAngle2) then
              ! add the contribution to the potential
              mass = gr_bhTreeArray(cpu, tr)%p(pos)
              if ( (gr_bhBndType(1) .EQ. GR_TREE_BND_PERIODIC) &
              .or. (gr_bhBndType(3) .EQ. GR_TREE_BND_PERIODIC) &
              .or. (gr_bhBndType(5) .EQ. GR_TREE_BND_PERIODIC)) then
                potential = potential + mass * gr_bhEwald(abs(dx), abs(dy), abs(dz))
              else
                potential = potential + mass*sqrt(dist2i)
              endif
              nodecnt = nodecnt + 1 ! distance counter
            else
              ! put all children on the stack
              mi(:) = multi(:)
              do i = 1,8
                mi(level+1) = i
                sp = sp + 1
                stack(:,sp) = mi(:)
              enddo
            endif
          endif
              
10        continue ! end 
          ! stack is empty - exiting
          if (sp == 0) exit
        enddo

        Phi(ii,jj,kk) = Phi(ii,jj,kk) + potential

      enddo
    enddo
  enddo
  return
  return
end

