!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhPotentialBlock
!!
!! NAME
!!
!!  gr_bhPotential
!!
!!
!! SYNOPSIS
!!
!!   gr_bhPotentialBlock(
!!          integer,intent(in) :: block,
!!          integer,intent(in) :: idensvar,
!!          integer,intent(in) :: ipotvar
!!          )
!!
!! DESCRIPTION
!!
!!   Computes the gravitational potential in a specific block.
!!
!! ARGUMENTS
!!
!!  block    - ID of a block where the potential is calculated
!!  idensvar - number of grid varible with density (recently not used)
!!  ipotvar  - number of grid varible with grav. potential
!!
!!***



subroutine gr_bhPotentialBlock(block, idensvar, ipotvar)

  use Grid_interface
  use gr_bhInterface, ONLY : gr_bhGetTreeSize,  gr_bhBlockRelationship
  use gr_bhLocalInterface, ONLY : gr_bhILContrib
  use gr_bhData, ONLY : GR_TREE_BND_PERIODIC, gr_bhBndType, &
    gr_bhLx, gr_bhLy, gr_bhLz, gr_bhTreeLevels, &
    gr_bhTreeBS, gr_bhLocCoords, gr_bhTreeMyPE, &
    gr_bhTreeLrefine, gr_bhTreeDiag2, gr_bhIlist, &
    gr_bhGravFac, gr_bhTreeDcount, gr_bhTreeLimAngle2, &
    gr_bhTreeNFirstLev, gr_bhTreeChild, gr_bhTreeNodetype, &
    gr_bhTreeFirstLevBlocks, gr_bhTreeParentTree
  use tree, ONLY : nodetype, child, mchild

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer, intent(IN) :: block, idensvar, ipotvar
  integer       :: i, j, k, cpu, tr, sp, ii, jj, kk, fac, l, level, pos, temp
  real          :: dcount(1:5), cellcnt, nodecnt
  integer       :: multi(1:gr_bhTreeLevels)
  real, POINTER, DIMENSION(:,:,:,:) :: solnData
  real, dimension(1:gr_bhTreeBS, 1:gr_bhTreeBS, 1:gr_bhTreeBS) :: Phi
  integer :: stack(2, MAXBLOCKS)
  real    :: disti, dist2i, dx, dy, dz, dvol, mass, x, y, z, potential, xi2
  integer, dimension(2,MDIM)   :: blkLimits,blkLimitsGC
  logical       :: gcell = .false.


  ! get information about the block
  call Grid_getBlkIndexLimits(block,blkLimits,blkLimitsGC) ! indexes necessary to write values into solnData
  call Grid_getBlkPtr(block,solnData,CENTER) ! pointer to the density and gp field
  
  ! reset distance counter
  do i = 1,5
    dcount(i) = 0.0
  enddo

  ! reset potential values to zeros
  do k = 1,gr_bhTreeBS
    do j = 1,gr_bhTreeBS
      do i = 1,gr_bhTreeBS
        Phi(i,j,k) = 0.0
      enddo
    enddo
  enddo

  ! put 1st level blocks on stack
  stack(1:2,1:gr_bhTreeNFirstLev) = gr_bhTreeFirstLevBlocks(1:2,1:gr_bhTreeNFirstLev)
  sp = gr_bhTreeNFirstLev

  ! walk the amr tree
  do

    ! pop a block from the stack
    tr = stack(1,sp)
    cpu  = stack(2,sp)
    sp = sp - 1

    ! check if interaction list can be used
    if (gr_bhIlist .ne. 0) then
      temp = gr_bhBlockRelationship(block, tr, cpu)
    else
      temp = -1
    endif
    if ((gr_bhIlist .ne. 0) .and. (temp .gt. 0) .and. (gr_bhTreeNodetype(tr,cpu) .eq. 1) &
    .and. (gr_bhTreeLrefine(tr, cpu) .eq. gr_bhTreeLrefine(block, gr_bhTreeMyPE))) then
      call gr_bhILContrib(block, tr, cpu, temp, Phi, cellcnt, nodecnt)
      dcount(4) = dcount(4) + cellcnt
      dcount(5) = dcount(5) + nodecnt

    else
      ! no interaction list
      if (gr_bhTreeNodetype(tr, cpu) .eq. 1) then ! 1 = LEAF block
        ! contribution of a LEAF block
        call gr_bhLeafContrib(block, tr, cpu, Phi, cellcnt, nodecnt)
        dcount(1) = dcount(1) + cellcnt
        dcount(2) = dcount(2) + nodecnt
      else
        ! contribution of a parent block
 
        ! determine the distance between this block and the parent block (tr, cpu)
        dx = abs(gr_bhLocCoords(gr_bhTreeBS+1, 1, block) - gr_bhTreeParentTree(2, tr, cpu)) &
        &  - 0.5*abs(gr_bhLocCoords(gr_bhTreeBS,1,block) - gr_bhLocCoords(1,1,block))
        dy = abs(gr_bhLocCoords(gr_bhTreeBS+1, 2, block) - gr_bhTreeParentTree(3, tr, cpu)) &
        &  - 0.5*abs(gr_bhLocCoords(gr_bhTreeBS,2,block) - gr_bhLocCoords(1,2,block))
        dz = abs(gr_bhLocCoords(gr_bhTreeBS+1, 3, block) - gr_bhTreeParentTree(4, tr, cpu)) &
        &  - 0.5*abs(gr_bhLocCoords(gr_bhTreeBS,3,block) - gr_bhLocCoords(1,3,block))
        dx = max(dx, 0.0)
        dy = max(dy, 0.0)
        dz = max(dz, 0.0)
        ! periodic boundaries
        if (gr_bhBndType(1) .EQ. GR_TREE_BND_PERIODIC) then
          dx = min(abs(dx), abs(dx+gr_bhLx), abs(dx-gr_bhLx))
        endif
        if (gr_bhBndType(3) .EQ. GR_TREE_BND_PERIODIC) then
          dy = min(abs(dy), abs(dy+gr_bhLy), abs(dy-gr_bhLy))
        endif
        if (gr_bhBndType(5) .EQ. GR_TREE_BND_PERIODIC) then
          dz = min(abs(dz), abs(dz+gr_bhLz), abs(dz-gr_bhLz))
        endif
        
        dist2i = 1.0/(dx*dx + dy*dy + dz*dz + 1.d-99)
        if ((gr_bhTreeDiag2(gr_bhTreeLrefine(tr,cpu))*dist2i < gr_bhTreeLimAngle2)) then ! add to potential or use children?
          call gr_bhParentContrib(block, tr, cpu, Phi)
          dcount(3) = dcount(3) + gr_bhTreeBS**3 ! "distance to parent block" counter increases just by number of cell in this block
        else
          ! put all children on the stack
          do ii = 1,mchild
            sp = sp + 1
            stack(1,sp) = gr_bhTreeChild(1,ii,tr,cpu)
            stack(2,sp) = gr_bhTreeChild(2,ii,tr,cpu)
          enddo

        endif ! add to potential or use children?
        
      endif ! LEAF block
    endif ! interaction list


    ! stack is empty - exiting
    if (sp == 0) exit
  enddo

  ! multiply potential with the gravity constant
  do k = 1,gr_bhTreeBS
    do j = 1,gr_bhTreeBS
      do i = 1,gr_bhTreeBS
         solnData(ipotvar, i-1+blkLimits(LOW,IAXIS), j-1+blkLimits(LOW,JAXIS), k-1+blkLimits(LOW,KAXIS)) = &
              -gr_bhGravFac*Phi(i,j,k)
      enddo
    enddo
  enddo
  ! release the block pointer
  call Grid_releaseBlkPtr(block,solnData, CENTER)
  gr_bhTreeDcount = gr_bhTreeDcount + dcount

  return
end subroutine gr_bhPotentialBlock

