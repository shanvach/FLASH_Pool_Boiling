!!****if* source/Grid/GridMain/Samrai/Grid_getListOfBlocks
!!
!! NAME
!!  Grid_getListOfBlocks
!!
!! SYNOPSIS
!!  #include "Flash.h"
!!
!!  Grid_getListOfBlocks(integer(IN)  :: blockType,
!!                       integer(OUT) :: listofBlocks(MAXBLOCKS), 
!!                       integer(OUT) :: count)
!!  
!! DESCRIPTION 
!!  Returns a list and the number of blocks of a specified type on the local processor
!!  
!!
!! ARGUMENTS
!!  blockType - Valid Values LEAF/PARENT/ANCESTOR
!!  listofBlocks - returned array holding the integer block number of all the blocks
!!                 on a local processor of type 'blockType'
!!  count - number of blocks returned in listofBlocks
!!
!!***


subroutine Grid_getListOfBlocks(blockType, listOfBlocks, count)

  use Grid_data, ONLY : gr_maxPatches
  
implicit none
#include "Flash.h"

  integer, intent(in) :: blockType
  integer,dimension(:),intent(out) :: listOfBlocks
  integer,intent(out) :: count

  !DEV: is there a notion of MAX_PATCHES per proc in SAMRAI?  
  !need to dimension level and patch, not sure how much
  integer, dimension(MAXBLOCKS) :: levelNum, patchNum
  
  call samrai_get_list_of_blocks(patchNum,levelNum, count)

  !Samrai keeps track of patches with a patch and level number
  !Flash just uses a single number, combine levelNum and patchNum to get
  !a local flash block number

  !DEV: verify starting index.  samrai starts counting at 0

  do i=1, count
     listOfBlocks(i) = levelNum(i) * gr_maxPatches + patchNum(i)
  end do

  return
end subroutine Grid_getListOfBlocks
