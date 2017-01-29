!!****if* source/Grid/GridMain/Samrai/Grid_getBlkSize
!!
!! NAME
!!  Grid_getBlkPhysicalSize
!!
!! SYNOPSIS
!!  #include "Flash.h"
!!
!!  Grid_getBlkPhysicalSize(integer(IN)  :: blockId,
!!                    integer(OUT) :: blockSize(MDIM))
!!  
!! DESCRIPTION 
!!  Gets the size of the block along each dimension by calculation
!!  the delta or each direction times the number of zones.  It is important
!!  to calculate the size of a block consistently because rounding errors
!!  can occur in some cases.  Use this function to get the size of the 
!!  block rather than calculating a new way
!!  
!! ARGUMENTS
!!  blockId -local block number
!!  blockSize - returned array of size MDIM holding the size of 
!!              each dimension of the block
!!
!!
!!***

subroutine Grid_getBlkPhysicalSize(blockId, blockSize)

  use Grid_data, ONLY : gr_maxPatches

  implicit none

#include "constants.h"

  integer,intent(in) :: blockId
  real,dimension(MDIM),intent(out) :: blockSize

  integer :: levelNum, patchNum
  
  !translate blockID to patchNum and levelNum
  levelNum = blockID / gr_maxPatches

  patchNum = mod(blockID, gr_maxPatches)

  call samrai_get_patch_phys_size(patchNum, levelNum, blockSize)

  return
end subroutine Grid_getBlkPhysicalSize
