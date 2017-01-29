!!****if* source/Grid/GridMain/Samrai/Grid_getBlkRefineLevel
!!
!! NAME
!!  Grid_getBlkRefineLevel
!!
!! SYNOPSIS
!!
!!
!!  Grid_getBlkRefineLevel(integer(IN)  :: blockId,
!!                      integer(OUT) :: refineLevel)
!!  
!! DESCRIPTION 
!!  Get the refinement level of a given block as denoted by blockId
!!
!! ARGUMENTS
!!  blockid - the local block number
!!  refineLevel - returned value, refinement level of block
!!
!!***

subroutine Grid_getBlkRefineLevel(blockId, refineLevel)
  
  use Grid_data, ONLY : gr_maxPatches

implicit none
  integer,intent(in) :: blockId
  integer,intent(out) :: refineLevel

  refineLevel = blockID / gr_maxPatches

  return
end subroutine Grid_getBlkRefineLevel














