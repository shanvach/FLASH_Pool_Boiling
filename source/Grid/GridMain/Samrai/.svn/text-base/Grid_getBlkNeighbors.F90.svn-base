!!****if* source/Grid/GridMain/Samrai/Grid_getBlkNeighbors
!!
!! NAME
!!  Grid_getBlkNeighbors
!!
!! SYNOPSIS
!!  #include "Flash.h"
!!
!!  Grid_getBlkNeighbors(integer(IN)  :: blockId,
!!                         integer(OUT) :: neigh(NDIM*4,NDIM+2))
!!                    
!! DESCRIPTION 
!!  Gets information about block's neighbors
!!  
!! ARGUMENTS 
!!
!!  blockId - the local blockId, as recognized by PM3
!!
!! DEV : The description may not be very accurate
!!  neigh - neigh(1,:) contains the local blockId of the block
!!        neigh(2,:) contains the processor it is on
!!        There is provision for two neighbors along each face. If there is
!!        only one neighbor along a face, then all entries for the second 
!!        neighbor along the face contain a NULL value
!!
!!***

subroutine Grid_getBlkNeighbors(blockId, neighbors)

  use Grid_data, ONLY : gr_maxPatches

  implicit none

#include "constants.h"

  integer, intent(in) :: blockId
  integer, dimension(2,MDIM*2),intent(out):: neighbors

  integer :: levelNum, patchNum
  
  !translate blockID to patchNum and levelNum
  levelNum = blockID / gr_maxPatches

  patchNum = mod(blockID, gr_maxPatches)

  !do we really need neighbors right now?  kda verify

  call samrai_get_patch_neighbors(patchNum, levelNum, neighbors)
  

  return
end subroutine Grid_getBlkNeighbors
