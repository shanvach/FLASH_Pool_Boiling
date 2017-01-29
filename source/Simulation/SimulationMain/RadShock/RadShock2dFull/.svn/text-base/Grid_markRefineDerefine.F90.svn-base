!!****if* source/Simulation/SimulationMain/RadShock/RadShock2dFull/Grid_markRefineDerefine
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
!! like, gr_meshMe or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

subroutine Grid_markRefineDerefine()

  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,&
                        gr_numRefineVars,gr_refine_var,gr_refineOnParticleCount,&
                        gr_enforceMaxRefinement, gr_maxRefine,&
                        gr_lrefineMaxRedDoByTime,&
                        gr_lrefineMaxRedDoByLogR,&
                        gr_lrefineCenterI,gr_lrefineCenterJ,gr_lrefineCenterK
  use tree, ONLY : newchild, refine, derefine, stay, nodetype, lnblocks, lrefine
  use Grid_interface, ONLY : Grid_fillGuardCells, Grid_markRefineSpecialized
  use Simulation_data

  implicit none

#include "constants.h"
#include "Flash.h"

  
  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,i,iref
  logical :: doEos=.true.
  integer,parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS
  logical,dimension(maskSize) :: gcMask
  real :: specs(7)

  if(gr_lrefineMaxRedDoByTime) then
     call gr_markDerefineByTime()
  end if

  ! that are implemented in this file need values in guardcells

  gcMask=.false.
  do i = 1,gr_numRefineVars
     iref = gr_refine_var(i)
     if (iref > 0) gcMask(iref) = .TRUE.
  end do

  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,doEos=.true.,&
       maskSize=maskSize, mask=gcMask, makeMaskConsistent=.true.,&
       selectBlockType=ACTIVE_BLKS)

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

#ifdef FLASH_GRID_PARAMESH2
  ! Make sure lrefine_min and lrefine_max are obeyed - KW
  if (gr_numRefineVars .LE. 0) then
     call gr_markRefineDerefine(-1, 0.0, 0.0, 0.0)
  end if
#endif

  if(gr_refineOnParticleCount)call gr_ptMarkRefineDerefine()

  if(gr_enforceMaxRefinement) call gr_enforceMaxRefine(3)

  if(gr_lrefineMaxRedDoByLogR) &
       call gr_unmarkRefineByLogRadius(gr_lrefineCenterI,&
       gr_lrefineCenterJ,gr_lrefineCenterK)

  if (sim_lrefmaxBase > 0) then

     ! Apply base lrefine_max:
     call gr_enforceMaxRefine(sim_lrefmaxBase)

     ! Apply specialized lrefine_max to Beryllium region:
     if(sim_lrefmaxBe > 0) then
        specs = 0
        specs(1) = sim_belrXMin
        specs(2) = sim_belrXMax
        specs(3) = sim_belrYMin
        specs(4) = sim_belrYMax
        specs(7) = 1.0     
        call Grid_markRefineSpecialized(RECTANGLE, 7, specs, sim_lrefmaxBe)
     end if

     ! Apply specialized lrefine_max to Polyimide region:
     if(sim_lrefmaxPoly > 0) then
        specs = 0
        specs(1) = sim_polylrXMin
        specs(2) = sim_polylrXMax
        specs(3) = sim_polylrYMin
        specs(4) = sim_polylrYMax
        specs(7) = 1.0     
        call Grid_markRefineSpecialized(RECTANGLE, 7, specs, sim_lrefmaxPoly)
     end if     
  end if

  return
end subroutine Grid_markRefineDerefine

