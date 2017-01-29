!!****if* source/Grid/GridMain/Samrai/Grid_getBlkType
!!
!! NAME
!!  Grid_getBlkType
!!
!! SYNOPSIS
!!
!!  Grid_getBlkType(integer(IN)  :: blockId,
!!                    integer(OUT) :: blockType)
!!  
!! DESCRIPTION 
!!  Get the type of Block, whether it is a leaf, a parent or
!!  higher up in the tree. If leaf, the solution advances on it
!!  always, if parent, sometimes else it doesn't
!!
!!***

subroutine Grid_getBlkType(blockId, blockType)
implicit none
  integer, intent(in)  :: blockId
  integer, intent(out) :: blockType
  return

  !DEV: kda, wait to see anshu's development of Grid_getBlkType
  !if this were only for leaf, parent nodes then this function would
  !be redundant, but if we can get the blocks on physical boundaries
  !then we might want it.

end subroutine Grid_getBlkType
