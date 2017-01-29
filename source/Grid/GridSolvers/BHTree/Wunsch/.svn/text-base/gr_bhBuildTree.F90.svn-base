!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhBuildTree
!!
!! NAME
!!
!!  gr_bhBuildTree
!!
!!
!! SYNOPSIS
!!
!!   gr_bhBuildTree()
!!
!! DESCRIPTION
!!
!!   Builds the tree for the Poisson solver. The first part builds trees 
!!   in all LEAF blocks. The second part distributes masses and mass 
!!   centres of all blocks to all CPUs (array gr_bhTreeParentTree). The 
!!   third part calculates masses and mass centres of all parent blocks
!!   (also in the array gr_bhTreeParentTree).
!!
!! ARGUMENTS
!!
!!
!!***



subroutine gr_bhBuildTree(idensvar)

  use Grid_interface, ONLY : Grid_getListOfBlocks
  use gr_bhData, ONLY : gr_bhTreeParentTree, gr_bhTreeMaxlnblocks, &
    gr_bhTreeNumProcs, gr_bhTreeArray, gr_bhTreeBS, gr_bhTreeLrefine, &
    gr_bhTreeChild, gr_bhTreeLnblocks, gr_bhTreeMyPE, gr_bhTreeNodetype
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use gr_bhInterface, ONLY : gr_bhBuildTreeBlock
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use tree, ONLY : lrefine, child, mchild, nodetype, lnblocks, maxblocks_tr
      
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(IN) :: idensvar

  integer :: blockID, blockCount, cpu, tr, lrefine_max, lrefine_min, lev, i, ierr
  integer :: tc, cc
  integer :: blockList(gr_bhTreeMaxlnblocks)
  real    :: m, xm, ym, zm
  integer :: ipr
  logical       :: gcell = .false.
  real,dimension(1:gr_bhTreeBS) :: xCoords, yCoords, zCoords

  call Timers_start("build_tree")
  call Timers_start("bt1: det_block_trees")
  
  call RuntimeParameters_get('lrefine_max',lrefine_max)
  call RuntimeParameters_get('lrefine_min',lrefine_min)
  
  call Grid_getListOfBlocks(LEAF,blockList,blockCount) ! get list of LEAF blocks

  ! reset gr_bhLocactiveBlocks and gr_bhTreeArray
  do tr = 1, gr_bhTreeMaxlnblocks
    do cpu = 0, gr_bhTreeNumProcs-1
      nullify(gr_bhTreeArray(cpu, tr)%p) ! allocate crashes if gr_bhTreeArray()%p != NULL
    enddo
  enddo

  ! fill gr_bhTreeParentTree array with -1
  do cpu = 0, gr_bhTreeNumProcs-1
    do i=1,gr_bhTreeMaxlnblocks
      gr_bhTreeParentTree(1,i,cpu) = -1
      gr_bhTreeParentTree(2,i,cpu) = -1
      gr_bhTreeParentTree(3,i,cpu) = -1
      gr_bhTreeParentTree(4,i,cpu) = -1
    enddo
  enddo

  ! call BuildTreeBlock for each LEAF block
  do blockID = 1, blockCount
    call gr_bhBuildTreeBlock(idensvar, blockList(blockID))
  enddo

  call Timers_stop("bt1: det_block_trees")
  call Timers_start("bt2: com_leaf_blocks")
  ! gr_bhTreeParentTree includes masses and MC positions of all LEAF blocks
  ! we need to communicate them to compute masses and MC positions of parent blocks 
  ! (children can be on different CPU)
  call gr_bhComParentTree(-1)
  call Timers_stop("bt2: com_leaf_blocks")

  call Timers_start("bt3: com_parent_blocks")
  do lev = lrefine_max-1, 1, -1
    do cpu = 0, gr_bhTreeNumProcs-1
      do tr = 1, gr_bhTreeLnblocks(cpu)
        if ((gr_bhTreeNodetype(tr, cpu) /= 1) .and. (gr_bhTreeLrefine(tr, cpu) .eq. lev)) then
          m = 0.0
          xm = 0.0
          ym = 0.0
          zm = 0.0
        
          do i = 1,mchild
            tc = gr_bhTreeChild(1,i,tr,cpu)
            cc = gr_bhTreeChild(2,i,tr,cpu)
            if ((gr_bhTreeParentTree(1, tc, cc) == -1) .or. (gr_bhTreeParentTree(2, tc, cc) == -1) &
     &     .or. (gr_bhTreeParentTree(3, tc, cc) == -1) .or. (gr_bhTreeParentTree(4, tc, cc) == -1)) then
              print *, "missing PT, myPE, tr, cpu = ", gr_bhTreeMyPE, tr, cpu, tc, cc
              print *, gr_bhTreeParentTree(1, tc, cc),gr_bhTreeParentTree(2, tc, cc), &
                       gr_bhTreeParentTree(3, tc, cc),gr_bhTreeParentTree(4, tc, cc) 
              call Driver_abortFlash("missing Parent Tree :(")
            endif
            m  =  m + gr_bhTreeParentTree(1, tc, cc)
            xm = xm + gr_bhTreeParentTree(2, tc, cc)*gr_bhTreeParentTree(1, tc, cc)
            ym = ym + gr_bhTreeParentTree(3, tc, cc)*gr_bhTreeParentTree(1, tc, cc)
            zm = zm + gr_bhTreeParentTree(4, tc, cc)*gr_bhTreeParentTree(1, tc, cc)
          enddo
          gr_bhTreeParentTree(1, tr, cpu) = m
          gr_bhTreeParentTree(2, tr, cpu) = xm/m
          gr_bhTreeParentTree(3, tr, cpu) = ym/m
          gr_bhTreeParentTree(4, tr, cpu) = zm/m
        endif
      enddo
    enddo
  enddo

  call Timers_stop("bt3: com_parent_blocks")
  
  call Timers_stop("build_tree")
  
  return
end subroutine gr_bhBuildTree

