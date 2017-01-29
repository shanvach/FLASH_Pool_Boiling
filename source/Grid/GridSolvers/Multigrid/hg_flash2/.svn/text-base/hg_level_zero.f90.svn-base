!******************************************************************************

! Routine:      hg_level_zero

! Description:  Zero a variable on a given level.

! Parameters:   ivar       Variable index for variable to zero
!               level      Refinement level of blocks to zero
!               LeafFlag   Flag controlling whether leaf blocks are zeroed.


subroutine hg_level_zero(level, ivar, LeafFlag)

!==============================================================================

use dBase, ONLY:  dBaseRefinementLevel, dBaseGetDataPtrSingleBlock, &
                  dBaseReleaseDataPtrSingleBlock, GC, dBasePropertyInteger, &
                  dBaseNodeType
use mg_common

implicit none

integer, intent(in) :: ivar, level, LeafFlag

integer             :: b, lnblocks
real, pointer       :: solnData(:,:,:,:)
logical             :: ZeroThisBlock

!==============================================================================

lnblocks = dBasePropertyInteger("LocalNumberOfBlocks")

do b = 1, lnblocks

  ZeroThisBlock = (dBaseRefinementLevel(b) == level)
  if (LeafFlag == MG_NODES_LEAF_ONLY) then
    ZeroThisBlock = (ZeroThisBlock .and. (dBaseNodeType(b) == 1))
  else if (LeafFlag == MG_NODES_PARENT_ONLY) then
    ZeroThisBlock = (ZeroThisBlock .and. (dBaseNodeType(b) /= 1))
  endif

  if (ZeroThisBlock) then

    solnData => dBaseGetDataPtrSingleBlock(b, GC)

    solnData(ivar,:,:,:) = 0.

    call dBaseReleaseDataPtrSingleBlock(b, solnData)

  endif

enddo

!==============================================================================

return
end subroutine hg_level_zero
