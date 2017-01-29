!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/ChooseLevel/gr_pfftValidateSelectedLevel
!!
!! NAME
!!
!!  gr_pfftValidateSelectedLevel
!!
!! SYNOPSIS
!!  
!!  gr_pfftValidateSelectedLevel(inOutLevel)
!!
!! DESCRIPTION
!!
!!  This subroutine ensures that the selected level can be computed 
!!  in ChooseLevel mode.  In ChooseLevel mode data is restricted / prolonged
!!  from LEAF block level to the chosen level.  As such, we must ensure that 
!!  the LEAF blocks are not resolved too finely relative to the solve level.
!!
!!  Example 1D numerical problem:
!!
!!  Blocks have NXB=8 and lrefine_max=6 and chosen level=2.
!!    At finest resolution the global grid contains 8 * 2**(6-1) = 256 cells.
!!    At level 2 the global grid contains 16 cells.
!!    If a LEAF block exists at level=6 its data must be restricted to level=2.
!!    This does not work in a straight forward way because:
!!       Data at level 6: 8 cells.
!!       Data at level 5: Restricted to 4 cells.
!!       Data at level 4: Restricted to 2 cells.
!!       Data at level 3: Restricted to 1 cells.
!!       Data at level 2: Restricted to 0 cells.
!!
!!***

subroutine gr_pfftValidateSelectedLevel(inOutLevel)

#include "constants.h"
#include "Flash.h"
#ifdef FLASH_GRID_PARAMESH
use tree, ONLY : lrefine_max
#endif
implicit none
integer, intent(INOUT) :: inOutLevel
integer :: deltaRefFactor, minCellsOnBlockSide
integer :: inLevel

#ifdef FLASH_GRID_PARAMESH
!Quick, dirty check is to look at lrefine_max, but we should eventually 
!cycle over all leaf blocks to obtain the max refinement in existance.  
!Replace all occurences of lrefine_max with this value.

inLevel = inOutLevel
deltaRefFactor = 2 ** (lrefine_max - inLevel)

minCellsOnBlockSide = NXB
if (NDIM > 1) minCellsOnBlockSide = min(minCellsOnBlockSide, NYB)
if (NDIM > 2) minCellsOnBlockSide = min(minCellsOnBlockSide, NZB)

if (deltaRefFactor > minCellsOnBlockSide) then
   inOutLevel = lrefine_max - &
        int(log(real(minCellsOnBlockSide)+tiny(1.0e0)) / log(2.0e0)) + 1
   print *, "Chosen level of:", inLevel, &
        "is too coarse.  Increasing level to:", inOutLevel
else
   inOutLevel = inLevel
end if
#endif

end subroutine gr_pfftValidateSelectedLevel
