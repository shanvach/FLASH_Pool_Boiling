!!***f* source/Grid/Grid_markBlkDerefine
!!
!! NAME
!!  Grid_markBlkDerefine
!!
!! SYNOPSIS
!!
!!  Grid_markBlkDerefine(integer : block,
!!                         logical : mark)
!!  
!! DESCRIPTION 
!!  Mark a block for refinement or derefinement
!!
!!
!! 
!! 
!! 
!! 
!!
!!***

subroutine Grid_markBlockDeRefine(block, mark)
  
  use Grid_data, ONLY : gr_maxPatches
  
implicit none
  integer, intent(IN) :: block
  logical, intent(IN) :: mark

  integer :: levelNum, patchNum
  
  !translate blockID to patchNum and levelNum
  levelNum = blockID / gr_maxPatches

  patchNum = mod(blockID, gr_maxPatches)

  call samrai_mark_patch_derefine()

end subroutine Grid_markBlockForRefineDerefine
