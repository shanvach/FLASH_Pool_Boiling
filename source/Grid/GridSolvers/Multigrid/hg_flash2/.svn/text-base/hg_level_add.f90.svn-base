!******************************************************************************

! Routine:      hg_level_add

! Description:  Add two mesh variables in all blocks on a given level, placing
!               the result in the first variable.  Boundary conditions are
!               left inconsistent for the first variable.

! Parameters:   ivar1      Variable index for first variable (receives result)
!               ivar2      Variable index for second variable (unchanged)
!               level      Refinement level of blocks to add
!               LeafFlag   Flag controlling whether leaf blocks are added.


subroutine hg_level_add(level, ivar1, ivar2, LeafFlag)

!==============================================================================

use dBase, ONLY:  dBaseRefinementLevel, dBaseGetDataPtrSingleBlock, &
                  dBaseReleaseDataPtrSingleBlock, GC, dBasePropertyInteger, &
                  dBaseNodeType
use mg_common

implicit none

integer, intent(in) :: ivar1, ivar2, level, LeafFlag

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

    solnData(ivar1,:,:,:) = solnData(ivar1,:,:,:) + solnData(ivar2,:,:,:)

    call dBaseReleaseDataPtrSingleBlock(b, solnData)

  endif

enddo

!==============================================================================

return
end subroutine hg_level_add
