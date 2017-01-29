!!****if* source/Grid/GridSolvers/Multigrid_forDiffuseAdvanceByHgADI/gr_hgSetMaxLevel
!!
!! NAME
!!  gr_hgSetMaxLevel
!!
!! SYNOPSIS
!!   gr_hgSetMaxLevel(integer, intent(in):: level)
!!
!! DESCRIPTION
!!  
!!  Reworks the nodetypes(:) array in paramesh such that blocks of
!!  lrefine = level are leaves, as well as leaves above level.  
!!  This then resets Paramesh to believe that this is the entire mesh
!!
!!  This work used to be done by gr_hgBndry; however this was called
!!  far too often to be efficient; the dataflow was also convoluted,
!!  It was moved here for these reasons.
!!
!! ARGUMENTS
!!  
!!  level - the level that is set as the max level in paramesh
!!
!! NOTES
!!
!!  We can define CAUTIOUS_PARAMESH_CALL to use the more expensive
!!  amr_morton_process().  We have run more tests using this call with 
!!  Multigrid.
!! 
!!***

#include  "Flash.h"
#include  "Multigrid.h"
#include  "constants.h"

subroutine gr_hgSetMaxLevel (level)

  use Grid_data, ONLY : gr_meshNumProcs, gr_meshMe
  use Grid_interface, ONLY : Grid_getNumProcs, Grid_getMyPE
  use tree, ONLY : nodetype, lrefine, lnblocks, newchild
  use gr_hgData, ONLY: gr_hgSaveNodetype, gr_hgSaveNewchild

  use Timers_interface, ONLY: Timers_start, Timers_stop
!!  use gr_hgInterface, ONLY: gr_hgBndry

  implicit none

  integer,intent(IN) :: level

  integer            :: lb, nprocs, myPE


  call Grid_getNumProcs(nprocs)
  call Grid_getMyPE(myPE)
  newchild(:) = .false.
  nodetype(:) = gr_hgsaveNodetype(:)
  do lb = 1, lnblocks
     if (lrefine(lb) > level)  nodetype(lb) = -1 !invalid
     if (lrefine(lb) == level) then 
        nodetype(lb) = LEAF
        newchild(lb) = .TRUE.
     endif
     if ((lrefine(lb) == level-1) .and. (nodetype(lb) /= LEAF)) & 
          nodetype(lb) = PARENT_BLK

  enddo

  call Timers_start("amr_morton_process")

#ifdef FLASH_GRID_PARAMESH2
  call get_tree_nodetypes (gr_meshNumProcs, gr_meshMe)
#else

#ifdef CAUTIOUS_PARAMESH_CALL
  call amr_morton_process()
#else
  call amr_get_new_nodetypes(nprocs, myPE, level)
#endif

#endif

  call Timers_stop("amr_morton_process")

end subroutine gr_hgSetMaxLevel
