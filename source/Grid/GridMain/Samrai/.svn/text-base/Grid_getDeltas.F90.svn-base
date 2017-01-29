!!****if* source/Grid/GridMain/Samrai/Grid_getDeltas
!!
!! NAME
!!  Grid_getDeltas
!!
!! SYNOPSIS
!!
!!  Grid_getDeltas(integer(IN) :: blockId,
!!                  real(OUT)   :: del(MDIM))
!!  
!! DESCRIPTION 
!!  
!!  Gets the dx/dy/dz for a given blockId on the Grid
!!  dx is the size of one zone in the x direction of a block
!!
!!  
!! ARGUMENTS 
!!
!!  blockId - local block number
!!  del - array of size MDIM returned holding the dx, dy, and dz values
!!
!!***

subroutine Grid_getDeltas(blockId, del)

  use Grid_data, ONLY : gr_maxPatches

  implicit none

#include "constants.h"
  integer, intent(IN)               :: blockId
  real, dimension(MDIM),intent(OUT) :: del

  integer :: levelNum, patchNum
  
  !translate blockID to patchNum and levelNum
  levelNum = blockID / gr_maxPatches

  patchNum = mod(blockID, gr_maxPatches)

  call samrai_get_deltas(patchNum, levelNum, del)

  return
end subroutine Grid_getDeltas

