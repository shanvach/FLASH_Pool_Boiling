!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhFinalize
!!
!! NAME
!!
!!  gr_bhFinalize
!!
!! 
!! SYNOPSIS
!!
!!  call gr_bhFinalize()
!!
!!
!! DESCRIPTION
!!
!!  Finalizes the tree Poisson solver.  Deallocates all arrays.
!!
!!***

subroutine gr_bhFinalize()

  use gr_bhData, ONLY: gr_bhTreeEwald, &
       gr_bhTreeArray,gr_bhLocSentTreeLevels,gr_bhLocRecvTreeLevels,gr_bhTreeCellSize, &
       gr_bhTreeDiag2,gr_bhTreeBCen,gr_bhLocCoords, &
       gr_bhTreeParentTree,gr_bhLocParentTree, &
       gr_bhTreeNodetype,gr_bhTreeLrefine,gr_bhTreeChild,gr_bhTreeLnblocks, &
       gr_bhIlist, gr_bhTreeNeigh,gr_bhTreeSurbox, &
       gr_bhTreeCn_ind,gr_bhTreeCc_disti,gr_bhTreeCc_count,gr_bhTreeCc_ind, &
       gr_bhTreeBlocklist, gr_bhTreeFirstLevBlocks

  implicit none

  if (allocated(gr_bhTreeEwald)) deallocate(gr_bhTreeEwald)

  deallocate(gr_bhTreeArray)
  deallocate(gr_bhLocSentTreeLevels)
  deallocate(gr_bhLocRecvTreeLevels)
  deallocate(gr_bhTreeCellSize)
  deallocate(gr_bhTreeDiag2)
  deallocate(gr_bhTreeBCen)
  deallocate(gr_bhLocCoords)

  deallocate(gr_bhTreeParentTree)
  deallocate(gr_bhLocParentTree)

  deallocate(gr_bhTreeNodetype)
  deallocate(gr_bhTreeLrefine)
  deallocate(gr_bhTreeChild)
  deallocate(gr_bhTreeLnblocks)
  deallocate(gr_bhTreeBlocklist)
  deallocate(gr_bhTreeFirstLevBlocks)

  ! INTERACTION LISTS
  if (gr_bhIlist .ne. 0) then
    deallocate(gr_bhTreeNeigh)
    deallocate(gr_bhTreeSurbox)

    deallocate(gr_bhTreeCn_ind)
    deallocate(gr_bhTreeCc_disti)
    deallocate(gr_bhTreeCc_count)
    deallocate(gr_bhTreeCc_ind)
  endif

  return
end subroutine gr_bhFinalize
