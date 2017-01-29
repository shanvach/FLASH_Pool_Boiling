!!****if* source/Grid/GridSolvers/Multigrid_forDiffuseAdvanceByHgADI/gr_hgLevelAddScalar
!!
!! NAME
!!  gr_hgLevelAddScalar
!!
!! SYNOPSIS
!!
!!  gr_hgLevelAddScalar(integer, intent(in) :: level
!!                      integer, intent(in) :: ivar
!!                      real, intent(in)    :: scalar
!!                      integer, intent(in) :: LeafFlag
!!
!! DESCRIPTION
!!
!!  Adds a scalar to a subset of the blocks
!!
!! ARGUMENTS
!!
!!  level     - The refinement of the leaf blocks to add
!!  ivar      - the first variable index
!!  scalar    - the scalar to add
!!  LeafFlag  - determines which blocks are filled
!!              MG_NODES_LEAF_ONLY  only add the leaf blocks
!!              MG_NODES_PARENT_ONLY  only add the parent blocks
!!              MG_NODES_ALL_NODES  add all blocks
!!
!! RESULT
!!  The variable is uniformly increased by scalar
!!
!!***

!!REORDER(5): unk

subroutine gr_hgLevelAddScalar(level, ivar, scalar, LeafFlag)

!================================================================

  use tree, ONLY : lrefine,lnblocks,nodetype
  use physicaldata, ONLY: unk

  implicit none

#include "Multigrid.h"
#include "constants.h"

  integer, intent(in) :: ivar, level, LeafFlag
  real, intent(in)    :: scalar

  integer             :: b
  logical             :: AddThisBlock

!==============================================================================

  
  do b = 1, lnblocks
     
     AddThisBlock = (lrefine(b) == level)
     if (LeafFlag == MG_NODES_LEAF_ONLY) then
        AddThisBlock = (AddThisBlock .and. (nodetype(b) == LEAF))
     else if (LeafFlag == MG_NODES_PARENT_ONLY) then
        AddThisBlock = (AddThisBlock .and. (nodetype(b) /= LEAF))
     endif
     
     if (AddThisBlock) then
        unk(ivar,:,:,:,b) = unk(ivar,:,:,:,b) + scalar
     
     endif
     
  enddo

  !========================================================================

  return
end subroutine gr_hgLevelAddScalar
