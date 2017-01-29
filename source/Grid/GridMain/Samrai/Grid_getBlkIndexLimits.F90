!!****if* source/Grid/GridMain/Samrai/Grid_getBlkIndexLimits
!!
!! NAME
!!  Grid_getBlkIndexLimits
!!
!! SYNOPSIS
!!  #include "Flash.h"
!!
!!  call Grid_getBlkIndexLimits(integer(IN)  :: blockId,
!!                      integer(OUT) :: blkLimits(HIGH,MDIM),
!!                      integer(OUT) :: blkLimitsGC(HIGH,MDIM))
!!  
!! DESCRIPTION 
!!  
!! Returns the indicies of the upper and lower bounds of a block.
!! blkLimitsGC holds the entire dimension of the block (including guardcells)
!! blkLimits returns the indicies of the interior of the block.
!!
!! The first dimension holds the lower and upper bounds of the block.
!! The second dimension is of size MDIM to hold the upper and lower indicies of each dimension
!!  
!! ARGUMENTS 
!!
!!  blockId :: the local blockID
!!  blkLimits   :: the blkLimits of indices of the block that contain actual data
!!  blkLimitsGC  :: the dimensioning of the array containing the block
!!
!!             
!! EXAMPLES
!! 
!!  Take a 2 dimensional block of size NXB x NYB, with ghost cells along
!!  each dimension being GX and GY.  This block could be stored in an array of 
!!  size (NXB+2*GX,NYB+2*GY). For this array the returned values of blkLimits and blkLimitsGC
!!  are as follows.
!!
!!  blkLimitsGC(LOW, IAXIS) = 1 !lower bound index of the first cell in the entire block in the xdir
!!  blkLimitsGC(LOW, JAXIS) = 1
!!  blkLimitsGC(LOW, KAXIS) = 1
!!
!!  blkLimitsGC(HIGH, IAXIS) = NXB+2*GX !upper bound index of the last cell in the entire block in the xdir
!!  blkLimitsGC(HIGH, JAXIS) = NYB +2*GY
!!  blkLimitsGC(HIGH, KAXIS) = 1 !because there are only 2 dimensions
!!
!!
!!  blkLimits(LOW, IAXIS) = GX+1 !lower bound index of the first interior cell of a block in the xdir
!!  blkLimits(LOW, JAXIS) = GY+1 
!!  blkLimits(LOW, KAXIS) = 1
!!
!!  blkLimits(HIGH, IAXIS) = NXB+GX !the upper bound index of the last interior cell of a block in the xdir
!!  blkLimits(HIGH, JAXIS) = NYB+GY
!!  blkLimits(HIGH, KAXIS) = 1 !because there are only 2 dimensions
!!
!!***

subroutine Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC)

  use Grid_data, ONLY : gr_maxPatches, gr_iguard, gr_jguard, gr_kguard

#include "constants.h"

  implicit none

  integer,intent(IN) :: blockId
  integer,dimension(2,MDIM), intent(out) :: blkLimits, blkLimitsGC

  integer, dimension(2,MDIM) :: rangeSamrai
  integer :: levelNum, patchNum

  patchNum = mod(blockId, gr_maxPatches)
  levelNum = blockId/gr_maxPatches

  !have samrai just return the blkLimits, we can calculate blkLimitsGC in f90
  call samrai_get_patch_index_limits(patchNum, levelNum, rangeSamrai)

  !calculate the blkLimitsGC!
  !in samrai first interior cell is 0
  blkLimits(LOW,IAXIS) = rangeSamrai(1, IAXIS) + gr_iguard + 1
  blkLimits(HIGH,IAXIS) = rangeSamrai(2, IAXIS) + gr_iguard + 1

  blkLimits(LOW,JAXIS) = rangeSamrai(1, JAXIS) + gr_jguard + 1
  blkLimits(HIGH,JAXIS) = rangeSamrai(2, JAXIS) + gr_jguard + 1

  blkLimits(LOW,KAXIS) = rangeSamrai(1, KAXIS) + gr_kguard + 1
  blkLimits(HIGH,KAXIS) = rangeSamrai(2, KAXIS) + gr_kguard + 1



  blkLimitsGC(LOW,IAXIS) = blkLimits(LOW,IAXIS) - gr_iguard
  blkLimitsGC(LOW,IAXIS) = blkLimits(LOW,IAXIS) + gr_iguard

  blkLimitsGC(LOW,JAXIS) = blkLimits(LOW,JAXIS) - gr_jguard
  blkLimitsGC(LOW,JAXIS) = blkLimits(LOW,JAXIS) + gr_jguard

  blkLimitsGC(LOW,KAXIS) = blkLimits(LOW,KAXIS) - gr_kguard
  blkLimitsGC(LOW,KAXIS) = blkLimits(LOW,KAXIS) + gr_kguard



end subroutine Grid_getBlkIndexLimits
