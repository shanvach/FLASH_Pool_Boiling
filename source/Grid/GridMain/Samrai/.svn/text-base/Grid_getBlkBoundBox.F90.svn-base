!!****if* source/Grid/GridMain/Samrai/Grid_getBlkBoundBox
!!
!! NAME
!!  Grid_getBlkBoundBox
!!
!! SYNOPSIS
!!  #include "Flash.h"
!! 
!!  Grid_getBlkBoundBox(integer(IN)  :: blockId,
!!                      real(OUT) :: boundBox(2, MDIM))
!!   
!! DESCRIPTION 
!!  Gets the bounding box of the block identified by blockId
!!  the boundbox is defined as the lower left and upper right
!!  corners of the block
!!
!! ARGUMENTS
!!  blockId -local block number
!!  boundBox - returned array holding the boundBox coordinates in
!!            each dimension
!!
!!***


subroutine Grid_getBlkBoundBox(blockId,boundBox)

  use Grid_data, ONLY : gr_maxPatches

implicit none
#include "constants.h"
  integer, intent(in) :: blockId
  real, dimension(2, MDIM), intent(out) :: boundBox

  integer :: levelNum, patchNum


  !translate blockID to patchNum and levelNum
  levelNum = blockID / gr_maxPatches

  patchNum = mod(blockID, gr_maxPatches)

  call samrai_get_patch_bound_box(patchNum, levelNum, boundBox)

  !just the interior right?

  return
end subroutine Grid_getBlkBoundBox
