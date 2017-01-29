!******************************************************************************

! Routine:      hg_level_smultiply

! Description:  Multiply a mesh variable by a scalar quantity in all blocks on
!               a given level.  Boundary conditions are left inconsistent.

! Parameters:   ivar       Variable index for mesh variable
!               scalar     Scalar value by which to multiply
!               level      Refinement level of blocks to multiply
!               LeafFlag   Flag controlling whether leaf blocks are multiplied


subroutine hg_level_smultiply(level, ivar, scalar, LeafFlag)

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
logical             :: MultiplyThisBlock

!==============================================================================

lnblocks = dBasePropertyInteger("LocalNumberOfBlocks")

do b = 1, lnblocks

  MultiplyThisBlock = (dBaseRefinementLevel(b) == level)
  if (LeafFlag == MG_NODES_LEAF_ONLY) then
    MultiplyThisBlock = (MultiplyThisBlock .and. (dBaseNodeType(b) == 1))
  else if (LeafFlag == MG_NODES_PARENT_ONLY) then
    MultiplyThisBlock = (MultiplyThisBlock .and. (dBaseNodeType(b) /= 1))
  endif

  if (MultiplyThisBlock) then

    solnData => dBaseGetDataPtrSingleBlock(b, GC)

    solnData(ivar,:,:,:) = solnData(ivar,:,:,:) * scalar

    call dBaseReleaseDataPtrSingleBlock(b, solnData)

  endif

enddo

!==============================================================================

return
end subroutine hg_level_smultiply
