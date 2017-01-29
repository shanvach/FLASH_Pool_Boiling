!!****if* source/Grid/GridMain/Samrai/Grid_getBlkCenterCoords
!!
!! NAME
!!  Grid_getBlkCenterCoords
!!
!! SYNOPSIS
!!
!!  Grid_getBlkCenterCoords(integer(IN) :: blockId,
!!                      real(OUT)   :: blockCenter(MDIM))
!!                      
!!  
!! DESCRIPTION 
!!   Gets the coordinates of the center of the block identified by
!!   blockId.  Returns the coordinates in an array blockCenter
!!
!!
!! ARGUMENTS
!!  blockId - local id number of the block. for UG always 1
!!  blockCenter - returned array of size MDIM holding the blockCenter coords
!!
!!
!!  
!!
!!***


subroutine Grid_getBlkCenterCoords(blockId, blockCenter)

  use Grid_data, ONLY : gr_maxPatches

  implicit none

#include "constants.h"

  integer,intent(in) :: blockId
  real,dimension(MDIM),intent(out) :: blockCenter

  integer :: levelNum, patchNum
  
  !translate blockID to patchNum and levelNum
  levelNum = blockID / gr_maxPatches

  patchNum = mod(blockID, gr_maxPatches)

  call samrai_get_patch_ctr_coords(patchNum, levelNum, blockCenter)

end subroutine Grid_getBlkCenterCoords
