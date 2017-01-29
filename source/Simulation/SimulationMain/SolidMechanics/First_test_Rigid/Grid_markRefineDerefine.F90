!!****if* source/Simulation/SimulationMain/SolidMechanics/First_test_Rigid/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  Grid_markRefineDerefine()
!!  
!! DESCRIPTION 
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual 
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub.
!!
!! ARGUMENTS
!! 
!! NOTES
!!
!! Every unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. For Grid unit these variables begin with "gr_"
!! like, gr_myPE or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

subroutine Grid_markRefineDerefine()

#include "Flash.h"

#ifdef FLASH_GRID_PARAMESH
  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,&
                        gr_numRefineVars,gr_refine_var
  use tree, ONLY : newchild, refine, derefine, stay,lrefine_max
  use Grid_interface, ONLY : Grid_markRefineSpecialized,Grid_fillGuardCells, &
                             Grid_getListOfBlocks,Grid_getBlkRefineLevel
  use Simulation_data

  use ImBound_Data, only : ib_BlockMarker_flag

  implicit none

#include "constants.h"

  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,i,iref,refineLevel


  logical :: gcMask(NUNK_VARS)


  integer :: blockCount
  integer :: blockList(MAXBLOCKS) 
 
  integer :: blockID,lb

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  call Grid_fillGuardCells(CENTER,ALLDIR)

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

  do l = 1,gr_numRefineVars
     iref = gr_refine_var(l)
     ref_cut = gr_refine_cutoff(l)
     deref_cut = gr_derefine_cutoff(l)
     ref_filter = gr_refine_filter(l)
     call gr_markRefineDerefine(iref,ref_cut,deref_cut,ref_filter)
  end do


  do lb=1,blockCount

     blockID = blockList(lb)

     call Grid_getBlkRefineLevel(blockID,refineLevel)

!     write(*,*) 'Block=',blockID,refineLevel,ib_BlockMarker_flag(blockID)

     if ( ib_BlockMarker_flag(blockID) ) then ! Tru = there are markers in the block
        refine(blockID)   = .TRUE.
        !stay(blockID)     = .TRUE.
        derefine(blockID) = .FALSE.
     else                                     ! No Markers -> Derefine
        derefine(blockID) = .TRUE.
        refine(blockID)   = .FALSE.
        !stay(blockID)     = .FALSE.
     endif
  enddo

#ifdef FLASH_GRID_PARAMESH2
  ! Make sure lrefine_min and lrefine_max are obeyed - KW
  if (gr_numRefineVars .LE. 0) then
     call gr_markRefineDerefine(-1, 0.0, 0.0, 0.0)
  end if
#endif
  return
#endif

end subroutine Grid_markRefineDerefine
