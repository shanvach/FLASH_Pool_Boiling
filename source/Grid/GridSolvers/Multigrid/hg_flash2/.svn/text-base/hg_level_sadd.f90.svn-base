!******************************************************************************

! Routine:      hg_level_sadd

! Description:  Add a scalar quantity to a mesh variable in all blocks on
!               a given level.  Boundary conditions are left inconsistent.

! Parameters:   ivar       Variable index for mesh variable
!               scalar     Scalar value to add
!               level      Refinement level of blocks to add
!               LeafFlag   Flag controlling whether leaf blocks are added


subroutine hg_level_sadd(level, ivar, scalar, LeafFlag)

!==============================================================================

use dBase, ONLY:  dBaseRefinementLevel, dBaseGetDataPtrSingleBlock, &
                  dBaseReleaseDataPtrSingleBlock, GC, dBasePropertyInteger, &
                  dBaseNodeType
use mg_common

implicit none

integer, intent(in) :: ivar, level, LeafFlag
real, intent(in)    :: scalar

integer             :: b, lnblocks
real, pointer       :: solnData(:,:,:,:)
logical             :: AddThisBlock

!==============================================================================

lnblocks = dBasePropertyInteger("LocalNumberOfBlocks")

do b = 1, lnblocks

  AddThisBlock = (dBaseRefinementLevel(b) == level)
  if (LeafFlag == MG_NODES_LEAF_ONLY) then
    AddThisBlock = (AddThisBlock .and. (dBaseNodeType(b) == 1))
  else if (LeafFlag == MG_NODES_PARENT_ONLY) then
    AddThisBlock = (AddThisBlock .and. (dBaseNodeType(b) /= 1))
  endif

  if (AddThisBlock) then

    solnData => dBaseGetDataPtrSingleBlock(b, GC)

    solnData(ivar,:,:,:) = solnData(ivar,:,:,:) + scalar

    call dBaseReleaseDataPtrSingleBlock(b, solnData)

  endif

enddo

!==============================================================================

return
end subroutine hg_level_sadd
