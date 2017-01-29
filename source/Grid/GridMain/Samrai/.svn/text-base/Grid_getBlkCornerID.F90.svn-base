!!****if* source/Grid/GridMain/Samrai/Grid_getBlkCornerID
!!
!! NAME
!!  Grid_getBlkCornerID
!!
!! SYNOPSIS
!!  #include "Flash.h"
!!
!!  call Grid_getBlkCornerID(integer(IN)  :: blockId,
!!                           integer(OUT) :: Index(MDIM)
!!                           integer(OUT) :: stride(MDIM)
!!  
!! DESCRIPTION 
!!  
!! Returns the global integer indices of the start of the interior zone
!! of the block
!! and the stride of indices along each dimensions.
!! For example, in a 1 dimensional UG case with 2 blocks and nxb=8
!! The index for block 1 = 1 and the index for block 2 = 9 
!!  
!! ARGUMENTS 
!!
!!  blockId :: the local blockID
!!  Index   :: global integer indices of start of  block 
!!  stride  :: stride for indices
!!             
!! EXAMPLES 
!!  
!!
!!***

subroutine Grid_getBlkCornerID(blockId, index, stride)

  use Grid_data, ONLY : gr_maxPatches

  implicit none

#include "constants.h"
#include "Flash.h"
  integer,intent(IN) :: blockId
  integer,dimension(MDIM), intent(OUT) :: index, stride

  integer :: levelNum, patchNum

  !translate blockID to patchNum and levelNum
  levelNum = blockID / gr_maxPatches

  patchNum = mod(blockID, gr_maxPatches)

  call samrai_get_patch_corner_id(patchNum, levelNum, index)

  !translate between fortran and C indexing
  index(IAXIS) = index(IAXIS) + 1
  index(JAXIS) = index(JAXIS) + 1
  index(KAXIS) = index(KAXIS) + 1


end subroutine Grid_getBlkCornerID

